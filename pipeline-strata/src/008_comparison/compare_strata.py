#!/usr/bin/env python3
"""Compare pQTL effects between strata, using the cross-strata meta as reference.

Reads the staged tree written by src/006_staging.sh:

    <staged>/<seqid>/<group>.tsv.gz      e.g. 21231-3/{ever,never,meta}.tsv.gz

METHOD
------
1. DISCOVERY, once, from the reference meta -- never from a stratum.
   Selecting pQTLs in one stratum and then comparing to the other produces
   winner's curse: variants chosen for a large effect regress toward the mean
   in the other arm, manufacturing differences out of noise.

   Selecting on an inverse-variance-weighted meta is exactly unbiased for the
   difference. With V the sampling variances and W = 1/V_a + 1/V_b,

       Cov(beta_meta, beta_a - beta_b) = (1/W) [ V_a/V_a - V_b/V_b ] = 0

   identically, for any variances, because IVW weights are proportional to
   1/V. METAL runs here with SCHEME STDERR, i.e. IVW, and fixed-effect IVW is
   associative, so the two-stage construction (within study across strata, then
   across studies) equals one-stage over all study x stratum estimates.

   Blind spot: a pQTL present in only one stratum has its effect roughly halved
   in the meta and may miss the threshold. Such variants are reported as a
   SECONDARY tier, explicitly flagged, because their difference statistics are
   selection-biased and must not be quoted as tests.

2. INDEPENDENCE by LD clumping against the same reference the heritability
   module uses, with the HLA region dropped (extended LD makes clumping there
   unreliable). Falls back to distance-based greedy pruning with no LD panel.

3. COMPARISON by a two-sample z test on the difference:

       delta = beta_a - beta_b,  SE = sqrt(SE_a^2 + SE_b^2),  z = delta / SE

   Valid because the strata are DISJOINT individuals, so the estimates are
   independent and variances add. This is Altman & Bland's interaction test.
   It is stronger than a pooled SNP x stratum term, which would assume common
   covariate effects and a common residual variance across strata; here each
   stratum is fitted entirely separately.

4. A GLOBAL test first, per-locus second. At ~4,000 per stratum the difference
   test has 1/sqrt(2) the precision of the main effect, so it detects
   near-absence in one arm rather than modulation. Asking whether the z
   distribution as a whole departs from N(0,1) is far better powered than any
   single locus, and control proteins give an empirical null.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REQUIRED = ["CHR", "POS", "ID", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def read_header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def run(cmd: str, log: Path | None = None) -> subprocess.CompletedProcess:
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if log is not None:
        with log.open("a") as fh:
            fh.write(f"$ {cmd}\n{res.stdout}{res.stderr}\n")
    return res


def build_assoc(meta: Path, out: Path, hla: str) -> None:
    """Flatten the reference meta to a PLINK --clump input.

    SNP ids are chr:pos, matching the LD reference built by 004 and the
    convention in bin/run_herit.sh. Duplicate positions keep the smallest P.
    """
    try:
        hc, ha, hb = hla.split(":")
    except ValueError:
        hc, ha, hb = "", "0", "0"
    awk = r"""
      BEGIN{FS=OFS="\t"}
      NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
      {
        chr=$(h["CHR"]); gsub(/^[Cc][Hh][Rr]/,"",chr); pos=$(h["POS"])
        if (chr=="" || pos=="") next
        if (chr==hc && pos+0>=ha+0 && pos+0<=hb+0) next
        p=$(h["P"]); if (p=="" || p=="NA") next
        snp=chr":"pos; p=p+0
        if (!(snp in seenP) || p < seenP[snp]) { seenP[snp]=p; row[snp]=snp OFS chr OFS pos OFS p }
      }
      END{ print "SNP","CHR","POS","P"; for(s in row) print row[s] }
    """
    # The program goes to a file rather than inline: repr()-ing a multi-line awk
    # script into a shell command turns its newlines into literal backslash-n.
    prog = out.with_suffix(".awk")
    prog.write_text(awk)
    subprocess.run(
        f"gzip -cd {str(meta)!r} | awk -v hc={hc!r} -v ha={ha!r} -v hb={hb!r} "
        f"-f {str(prog)!r} > {str(out)!r}",
        shell=True, check=True)


def clump(assoc: Path, workdir: Path, ld_ref: str, p1: float, r2: float,
          kb: int, log: Path) -> list[str]:
    """Independent index SNPs (chr:pos). LD clumping, else distance pruning."""
    if ld_ref and Path(f"{ld_ref}.bed").exists() and shutil.which("plink"):
        out = workdir / "clumped"
        res = run(
            f"plink --bfile {ld_ref!r} --clump {str(assoc)!r} "
            f"--clump-snp-field SNP --clump-field P "
            f"--clump-p1 {p1} --clump-p2 {p1} --clump-r2 {r2} --clump-kb {kb} "
            f"--out {str(out)!r}", log)
        clumped = Path(f"{out}.clumped")
        if clumped.exists():
            idx = [ln.split()[2] for ln in clumped.read_text().splitlines()[1:]
                   if ln.strip()]
            return [s for s in idx if s]
        with log.open("a") as fh:
            fh.write(f"plink --clump produced no output (rc={res.returncode})\n")
        return []

    # Distance-based greedy pruning: strongest first, keep if no kept index
    # within +/- kb on the same chromosome.
    df = pd.read_csv(assoc, sep="\t")
    df = df[df["P"] <= p1].sort_values("P")
    kept: list[tuple[str, int]] = []
    out_ids: list[str] = []
    window = kb * 1000
    for snp, chrom, pos in zip(df["SNP"], df["CHR"].astype(str), df["POS"].astype(int)):
        if any(c == chrom and abs(pos - p) <= window for c, p in kept):
            continue
        kept.append((chrom, pos))
        out_ids.append(snp)
    return out_ids


def extract_snps(path: Path, snp_ids: set[str], workdir: Path, tag: str) -> pd.DataFrame:
    """Pull just the index variants out of an 11M-row sumstats file.

    Filtering with awk keeps this to one streaming pass per file; reading the
    whole table into pandas for a few hundred rows would dominate the runtime.
    """
    if not snp_ids:
        return pd.DataFrame(columns=REQUIRED + ["SNP"])
    ids_file = workdir / f"ids_{tag}.txt"
    ids_file.write_text("\n".join(sorted(snp_ids)) + "\n")
    out = workdir / f"hit_{tag}.tsv"
    awk = r"""
      BEGIN{FS=OFS="\t"}
      NR==FNR { want[$1]; next }
      FNR==1  { for(i=1;i<=NF;i++) h[$i]=i
                print "SNP","CHR","POS","ID","EA","OA","EAF","BETA","SE","P","N","HET_ISQ"; next }
      {
        chr=$(h["CHR"]); gsub(/^[Cc][Hh][Rr]/,"",chr); pos=$(h["POS"])
        snp=chr":"pos
        if (!(snp in want)) next
        iq = (("HET_ISQ" in h) ? $(h["HET_ISQ"]) : "NA")
        print snp, chr, pos, $(h["ID"]), $(h["EA"]), $(h["OA"]), $(h["EAF"]),
              $(h["BETA"]), $(h["SE"]), $(h["P"]), $(h["N"]), iq
      }
    """
    prog = workdir / f"extract_{tag}.awk"
    prog.write_text(awk)
    subprocess.run(
        f"gzip -cd {str(path)!r} | awk -f {str(prog)!r} {str(ids_file)!r} - > {str(out)!r}",
        shell=True, check=True)
    df = pd.read_csv(out, sep="\t", dtype={"SNP": str, "CHR": str, "EA": str, "OA": str})
    # A position can appear more than once (multi-allelic); keep the strongest.
    if not df.empty:
        df = df.sort_values("P").drop_duplicates("SNP", keep="first")
    return df


def harmonise(ref: pd.DataFrame, alt: pd.DataFrame, suffix: str) -> pd.DataFrame:
    """Align alt to ref's effect allele, flipping BETA and EAF where needed.

    ever and never are separate METAL runs, so the same variant can come back
    with EA and OA swapped. Differencing BETA without aligning would flip signs
    at random and fabricate differences. Alleles are upper-cased first because
    METAL reports them in lower case.
    """
    a = alt.copy()
    for c in ("EA", "OA"):
        a[c] = a[c].astype(str).str.upper()
    m = ref[["SNP", "EA", "OA"]].copy()
    for c in ("EA", "OA"):
        m[c] = m[c].astype(str).str.upper()
    m = m.rename(columns={"EA": "EA_ref", "OA": "OA_ref"})

    a = a.merge(m, on="SNP", how="inner")
    same = (a["EA"] == a["EA_ref"]) & (a["OA"] == a["OA_ref"])
    flip = (a["EA"] == a["OA_ref"]) & (a["OA"] == a["EA_ref"])
    a["ALIGNED"] = same | flip

    a.loc[flip, "BETA"] = -a.loc[flip, "BETA"]
    a.loc[flip, "EAF"] = 1.0 - a.loc[flip, "EAF"]

    drop = ["EA", "OA", "EA_ref", "OA_ref", "ALIGNED", "CHR", "POS"]
    a = a[a["ALIGNED"]].drop(columns=[c for c in drop if c in a.columns])
    return a.rename(columns={c: f"{c}{suffix}" for c in
                             ("BETA", "SE", "P", "N", "EAF", "HET_ISQ", "ID")})


def annotate_cis(df: pd.DataFrame, protein: str, gene_map: pd.DataFrame | None,
                 window_kb: int) -> pd.DataFrame:
    """Label each index variant cis or trans relative to the encoding gene."""
    df = df.copy()
    df["CIS"] = "NA"
    df["GENE"] = "NA"
    if gene_map is None:
        return df
    row = gene_map[gene_map["protein_id"].astype(str) == str(protein)]
    if row.empty:
        return df
    r = row.iloc[0]
    gene_chr = str(r["chr"]).replace("chr", "")
    try:
        start, end = int(r["gene_start"]), int(r["gene_end"])
    except (TypeError, ValueError):
        return df
    w = window_kb * 1000
    near = (df["CHR"].astype(str) == gene_chr) & \
           (df["POS"].astype(int) >= start - w) & (df["POS"].astype(int) <= end + w)
    df["CIS"] = np.where(near, "cis", "trans")
    df["GENE"] = str(r.get("gene_symbol", "") or "NA")
    return df


def compare_pair(merged: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Two-sample z test on the between-stratum difference."""
    d = merged.copy()
    d["DELTA"] = d[f"BETA_{a}"] - d[f"BETA_{b}"]
    d["SE_DELTA"] = np.sqrt(d[f"SE_{a}"] ** 2 + d[f"SE_{b}"] ** 2)
    d["Z_DIFF"] = np.where(d["SE_DELTA"] > 0, d["DELTA"] / d["SE_DELTA"], np.nan)
    d["P_DIFF"] = 2.0 * stats.norm.sf(np.abs(d["Z_DIFF"]))
    # Cochran's Q over two groups is exactly z^2; reported because readers
    # expect it, not because it adds information.
    d["Q"] = d["Z_DIFF"] ** 2
    d["I2"] = np.clip((d["Q"] - 1) / d["Q"].replace(0, np.nan) * 100, 0, 100)
    d["SAME_DIRECTION"] = np.sign(d[f"BETA_{a}"]) == np.sign(d[f"BETA_{b}"])
    return d


def qc_flags(d: pd.DataFrame, a: str, b: str, max_isq: float,
             max_eaf_diff: float) -> pd.DataFrame:
    """Gates that must pass before a difference is believable.

    EAF_DIFF: the strata are the same populations, so allele frequencies should
    agree. A discrepancy means a strand or allele problem, which would surface
    as a spurious effect difference.

    HET_ISQ: high across-study heterogeneity WITHIN a stratum means that
    stratum's estimate is unstable, and a difference against the other may be
    one discordant cohort rather than the exposure.
    """
    d = d.copy()
    d["EAF_DIFF"] = (d[f"EAF_{a}"] - d[f"EAF_{b}"]).abs()
    isq = pd.concat([pd.to_numeric(d.get(f"HET_ISQ_{a}"), errors="coerce"),
                     pd.to_numeric(d.get(f"HET_ISQ_{b}"), errors="coerce")], axis=1)
    d["MAX_HET_ISQ"] = isq.max(axis=1)
    d["QC_EAF"] = d["EAF_DIFF"] <= max_eaf_diff
    d["QC_HET"] = d["MAX_HET_ISQ"].isna() | (d["MAX_HET_ISQ"] <= max_isq)
    d["QC_PASS"] = d["QC_EAF"] & d["QC_HET"]
    return d


def global_test(z: np.ndarray, label: str) -> dict:
    """Is the z distribution as a whole heavier than N(0,1)?

    Better powered than any single locus: a systematic between-stratum
    difference shifts the variance of z even when no individual variant clears
    a multiple-testing threshold.
    """
    z = z[np.isfinite(z)]
    if z.size < 3:
        return {"set": label, "n": int(z.size), "mean_z": np.nan, "var_z": np.nan,
                "ks_p": np.nan, "lambda": np.nan, "n_p_lt_05": 0, "frac_p_lt_05": np.nan}
    p = 2.0 * stats.norm.sf(np.abs(z))
    return {
        "set": label,
        "n": int(z.size),
        "mean_z": float(np.mean(z)),
        "var_z": float(np.var(z, ddof=1)),
        # Under H0 var(z)=1; >1 indicates real between-stratum heterogeneity.
        "ks_p": float(stats.kstest(z, "norm").pvalue),
        "lambda": float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1)),
        "n_p_lt_05": int((p < 0.05).sum()),
        "frac_p_lt_05": float((p < 0.05).mean()),
    }


def make_plots(loci: pd.DataFrame, outdir: Path, a: str, b: str) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    ok = loci[loci["QC_PASS"]].copy()
    if ok.empty:
        return False

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # 1. effect in one stratum against the other, with the line of equality
    ax = axes[0]
    for tier, mark in (("primary", "o"), ("secondary", "^")):
        s = ok[ok["TIER"] == tier]
        if not s.empty:
            ax.errorbar(s[f"BETA_{a}"], s[f"BETA_{b}"],
                        xerr=s[f"SE_{a}"], yerr=s[f"SE_{b}"],
                        fmt=mark, ms=4, lw=0.6, alpha=0.7, label=tier)
    lim = np.nanmax(np.abs(np.r_[ok[f"BETA_{a}"], ok[f"BETA_{b}"]])) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel(f"beta ({a})"); ax.set_ylabel(f"beta ({b})")
    ax.set_title("Effect by stratum"); ax.legend(fontsize=8)

    # 2. QQ of the difference p-values, controls overlaid as the empirical null
    ax = axes[1]
    for subset, colour, name in (
        (ok[ok["CONTROL"]], "tab:grey", "control proteins"),
        (ok[~ok["CONTROL"]], "tab:blue", "targeted proteins"),
    ):
        p = subset["P_DIFF"].dropna().values
        if p.size == 0:
            continue
        obs = -np.log10(np.sort(p))
        exp = -np.log10((np.arange(1, p.size + 1) - 0.5) / p.size)
        ax.plot(exp, obs, "o", ms=4, color=colour, label=f"{name} (n={p.size})")
    m = ax.get_xlim()[1]
    ax.plot([0, m], [0, m], "k--", lw=0.8)
    ax.set_xlabel("expected -log10(p)"); ax.set_ylabel("observed -log10(p)")
    ax.set_title("Difference test calibration"); ax.legend(fontsize=8)

    # 3. z distribution against the standard normal
    ax = axes[2]
    z = ok["Z_DIFF"].dropna().values
    if z.size:
        ax.hist(z, bins=max(10, min(40, z.size // 2)), density=True,
                alpha=0.6, color="tab:blue")
        xs = np.linspace(-4, 4, 200)
        ax.plot(xs, stats.norm.pdf(xs), "k--", lw=1.2, label="N(0,1)")
        ax.set_xlabel("z (difference)"); ax.set_ylabel("density")
        ax.set_title(f"var(z) = {np.var(z, ddof=1):.2f}"); ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(outdir / "comparison_plots.png", dpi=150)
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--staged", required=True, help="staged tree: <seqid>/<group>.tsv.gz")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--strata", default="ever,never", help="the two strata to compare")
    ap.add_argument("--reference", default="meta", help="group used for DISCOVERY only")
    ap.add_argument("--p-discovery", type=float, default=5e-8)
    ap.add_argument("--ld-ref", default=os.environ.get("GWAS_LD_REF", ""),
                    help="PLINK bfile prefix; empty => distance-based pruning")
    ap.add_argument("--clump-r2", type=float, default=0.1)
    ap.add_argument("--clump-kb", type=int, default=1000)
    ap.add_argument("--hla", default="6:25500000:34000000")
    ap.add_argument("--gene-map", default="", help="protein_gene_map.tsv for cis/trans")
    ap.add_argument("--cis-window-kb", type=int, default=1000)
    ap.add_argument("--annot", default="",
                    help="<panel>.annot.tsv; rows with source=random_control "
                         "form the empirical null")
    ap.add_argument("--max-het-isq", type=float, default=75.0)
    ap.add_argument("--max-eaf-diff", type=float, default=0.05)
    ap.add_argument("--proteins", default="", help="comma-separated subset")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    strata = [s.strip() for s in args.strata.split(",") if s.strip()]
    if len(strata) != 2:
        sys.exit(f"ERROR: --strata needs exactly two groups, got {strata}")
    a, b = strata
    ref = args.reference.strip()

    staged = Path(args.staged)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    log = outdir / "comparison.log"
    log.write_text("")

    gene_map = None
    if args.gene_map and Path(args.gene_map).exists():
        gene_map = pd.read_csv(args.gene_map, sep="\t", dtype=str)

    controls: set[str] = set()
    if args.annot and Path(args.annot).exists():
        ann = pd.read_csv(args.annot, sep="\t", dtype=str)
        controls = set(ann.loc[ann["source"] == "random_control", "seqid"])

    wanted = {p.strip() for p in args.proteins.split(",") if p.strip()}
    proteins = sorted(d.name for d in staged.iterdir()
                      if d.is_dir() and (not wanted or d.name in wanted))
    if not proteins:
        sys.exit(f"ERROR: no protein directories under {staged}")

    print(f"Comparing '{a}' vs '{b}', discovery from '{ref}' at P < {args.p_discovery:g}")
    print(f"{len(proteins)} protein(s); {len(controls)} flagged as controls\n")

    all_loci: list[pd.DataFrame] = []
    skipped: list[tuple[str, str]] = []

    for protein in proteins:
        pdir = staged / protein
        files = {g: pdir / f"{g}.tsv.gz" for g in (a, b, ref)}
        missing = [g for g, f in files.items() if not f.exists()]
        if missing:
            skipped.append((protein, f"missing {','.join(missing)}"))
            continue

        for g, f in files.items():
            hdr = read_header(f)
            absent = [c for c in REQUIRED if c not in hdr]
            if absent:
                skipped.append((protein, f"{g} lacks {','.join(absent)}"))
                break
        else:
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                assoc = work / "assoc.tsv"
                build_assoc(files[ref], assoc, args.hla)

                index = clump(assoc, work, args.ld_ref, args.p_discovery,
                              args.clump_r2, args.clump_kb, log)

                # SECONDARY tier: significant in a stratum but not in the
                # reference. Selection-biased by construction, so carried
                # separately and never reported as a test.
                sec_index: list[str] = []
                for g in (a, b):
                    sassoc = work / f"assoc_{g}.tsv"
                    build_assoc(files[g], sassoc, args.hla)
                    sec_index += clump(sassoc, work, args.ld_ref, args.p_discovery,
                                       args.clump_r2, args.clump_kb, log)
                sec_only = sorted(set(sec_index) - set(index))

                tiers = {s: "primary" for s in index}
                tiers.update({s: "secondary" for s in sec_only})
                if not tiers:
                    skipped.append((protein, "no variant reached the threshold"))
                    continue

                ids = set(tiers)
                dref = extract_snps(files[ref], ids, work, "ref")
                if dref.empty:
                    skipped.append((protein, "index variants absent from the reference"))
                    continue
                da = harmonise(dref, extract_snps(files[a], ids, work, a), f"_{a}")
                db = harmonise(dref, extract_snps(files[b], ids, work, b), f"_{b}")

                base = dref.rename(columns={c: f"{c}_{ref}" for c in
                                            ("BETA", "SE", "P", "N", "EAF", "HET_ISQ", "ID")})
                m = base.merge(da, on="SNP", how="inner").merge(db, on="SNP", how="inner")
                if m.empty:
                    skipped.append((protein, "no variant aligned across all three files"))
                    continue

                m = compare_pair(m, a, b)
                m = qc_flags(m, a, b, args.max_het_isq, args.max_eaf_diff)
                m = annotate_cis(m, protein, gene_map, args.cis_window_kb)
                m.insert(0, "PROTEIN", protein)
                m["TIER"] = m["SNP"].map(tiers)
                m["CONTROL"] = protein in controls
                all_loci.append(m)

                n_pri = int((m["TIER"] == "primary").sum())
                print(f"  {protein:<12} {n_pri:>4} primary, "
                      f"{len(m) - n_pri:>3} secondary index variant(s)")

    if not all_loci:
        for p, why in skipped:
            print(f"  {p:<12} skipped: {why}")
        sys.exit("\nNo protein yielded index variants; nothing to compare.")

    loci = pd.concat(all_loci, ignore_index=True)

    # Multiple testing over PRIMARY, QC-passing, non-control variants only:
    # those are the actual hypothesis tests.
    testable = loci["QC_PASS"] & (loci["TIER"] == "primary") & (~loci["CONTROL"])
    n_tests = int(testable.sum())
    loci["P_DIFF_BONF"] = np.where(testable,
                                   np.minimum(loci["P_DIFF"] * max(n_tests, 1), 1.0),
                                   np.nan)
    loci["SIGNIF_DIFF"] = loci["P_DIFF_BONF"] < 0.05

    loci = loci.sort_values(["PROTEIN", "P_DIFF"])

    front = ["PROTEIN", "SNP", "CHR", "POS", "EA", "OA", "GENE", "CIS",
             "TIER", "CONTROL"]
    stats_cols = ["DELTA", "SE_DELTA", "Z_DIFF", "P_DIFF", "P_DIFF_BONF",
                  "SIGNIF_DIFF", "SAME_DIRECTION", "Q", "I2"]
    per_group = [f"{m}_{g}" for g in (a, b, ref)
                 for m in ("BETA", "SE", "P", "N", "EAF", "HET_ISQ")]
    qc = ["EAF_DIFF", "MAX_HET_ISQ", "QC_EAF", "QC_HET", "QC_PASS"]
    order = [c for c in front + stats_cols + per_group + qc if c in loci.columns]
    order += [c for c in loci.columns if c not in order]
    loci = loci[order]

    loci.to_csv(outdir / "pqtl_comparison.tsv", sep="\t", index=False, na_rep="NA")

    # Global tests: the primary question, and the control set as empirical null.
    rows = [
        global_test(loci.loc[testable, "Z_DIFF"].values, "primary_targeted"),
        global_test(loci.loc[loci["QC_PASS"] & loci["CONTROL"], "Z_DIFF"].values,
                    "primary_controls"),
        global_test(loci.loc[loci["QC_PASS"] & (loci["TIER"] == "secondary"), "Z_DIFF"].values,
                    "secondary_selection_biased"),
    ]
    gt = pd.DataFrame(rows)
    gt.to_csv(outdir / "global_test.tsv", sep="\t", index=False, na_rep="NA")

    plotted = make_plots(loci, outdir, a, b)

    # ---- report ----------------------------------------------------------
    print(f"\n{'='*70}\nGLOBAL TEST -- is the z distribution heavier than N(0,1)?")
    print(f"{'set':<28}{'n':>5}{'mean z':>9}{'var z':>8}{'KS p':>10}{'p<.05':>8}")
    for r in rows:
        if r["n"] == 0:
            continue
        print(f"{r['set']:<28}{r['n']:>5}{r['mean_z']:>9.3f}{r['var_z']:>8.2f}"
              f"{r['ks_p']:>10.3g}{r['frac_p_lt_05']:>8.2f}")
    print("Under the null var(z) = 1. Values above 1 in the targeted set but not")
    print("the controls indicate genuine between-stratum heterogeneity.")

    n_fail = int((~loci["QC_PASS"]).sum())
    print(f"\n{'='*70}\nPER-LOCUS")
    print(f"  index variants            : {len(loci)}")
    print(f"  primary / secondary       : {int((loci['TIER']=='primary').sum())}"
          f" / {int((loci['TIER']=='secondary').sum())}")
    print(f"  failed QC (EAF or I2)     : {n_fail}")
    print(f"  hypothesis tests          : {n_tests}")
    hits = loci[loci["SIGNIF_DIFF"].fillna(False)]
    print(f"  Bonferroni-significant    : {len(hits)}")
    if not hits.empty:
        cols = ["PROTEIN", "SNP", "GENE", "CIS", f"BETA_{a}", f"SE_{a}",
                f"BETA_{b}", f"SE_{b}", "DELTA", "P_DIFF", "P_DIFF_BONF"]
        print("\n" + hits[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    if skipped:
        print(f"\nSkipped {len(skipped)} protein(s):")
        for p, why in skipped:
            print(f"  {p:<12} {why}")

    print(f"\nWrote {outdir/'pqtl_comparison.tsv'}")
    print(f"      {outdir/'global_test.tsv'}")
    if plotted:
        print(f"      {outdir/'comparison_plots.png'}")
    else:
        print("      (no plots: matplotlib unavailable)")


if __name__ == "__main__":
    main()
