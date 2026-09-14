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

# Decompression is the single largest cost. pigz parallelises it; fall back to
# gzip, which is single-threaded.
DECOMP = "pigz -dc" if shutil.which("pigz") else "gzip -cd"


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


def scan_sumstats(path: Path, workdir: Path, tag: str, p_max: float,
                  ids: set[str] | None, hla: str) -> pd.DataFrame:
    """ONE streaming pass, returning every row that is either below p_max OR at
    one of `ids`.

    The old code read each file twice: once to build a PLINK --clump input and
    again to pull the index variants out. Six full decompressions per protein
    (3 files x 2) is 3.5 TB for batch1. A single pass cannot be done with a
    p-filter alone -- an index variant chosen from the reference may be null in
    a stratum, so it has to be fetched by ID -- hence the OR.

    Returned rows carry the SNP key (chr:pos), so the caller can both clump and
    extract from this one frame.
    """
    try:
        hc, ha, hb = hla.split(":")
    except ValueError:
        hc, ha, hb = "", "0", "0"

    # The NR==FNR idiom below distinguishes the id file from the sumstats by
    # record count. An EMPTY id file contributes no records, so the sumstats
    # header would satisfy NR==FNR, be swallowed as an id, and leave the column
    # map empty. Always write a sentinel that cannot collide with a chr:pos key.
    ids_file = workdir / f"ids_{tag}.txt"
    ids_file.write_text("\n".join(sorted(ids) if ids else ["__NO_IDS__"]) + "\n")

    awk = r"""
      BEGIN{FS=OFS="\t"}
      NR==FNR { want[$1]; next }
      FNR==1  { for(i=1;i<=NF;i++) h[$i]=i
                print "SNP","CHR","POS","ID","EA","OA","EAF","BETA","SE","P","N","HET_ISQ"; next }
      {
        chr=$(h["CHR"]); gsub(/^[Cc][Hh][Rr]/,"",chr); pos=$(h["POS"])
        if (chr=="" || pos=="") next
        if (chr==hc && pos+0>=ha+0 && pos+0<=hb+0) next
        p=$(h["P"]); if (p=="" || p=="NA") next
        snp=chr":"pos
        if (!((snp in want) || (p+0 < pmax))) next
        iq = (("HET_ISQ" in h) ? $(h["HET_ISQ"]) : "NA")
        print snp, chr, pos, $(h["ID"]), $(h["EA"]), $(h["OA"]), $(h["EAF"]),
              $(h["BETA"]), $(h["SE"]), p, $(h["N"]), iq
      }
    """
    prog = workdir / f"scan_{tag}.awk"
    prog.write_text(awk)
    out = workdir / f"scan_{tag}.tsv"
    subprocess.run(
        f"{DECOMP} {str(path)!r} | awk -v hc={hc!r} -v ha={ha!r} -v hb={hb!r} "
        f"-v pmax={p_max:g} -f {str(prog)!r} {str(ids_file)!r} - > {str(out)!r}",
        shell=True, check=True)

    df = pd.read_csv(out, sep="\t", dtype={"SNP": str, "CHR": str, "EA": str, "OA": str})
    if not df.empty:
        # A position can appear more than once (multi-allelic); keep the strongest.
        df = df.sort_values("P").drop_duplicates("SNP", keep="first")
    return df


def clump_frame(df: pd.DataFrame, workdir: Path, tag: str, ld_ref: str,
                p1: float, r2: float, kb: int, log: Path,
                plink_mb: int = 2000) -> list[str]:
    """Independent index SNPs from an already-filtered frame.

    The frame passed here holds only variants below the threshold -- typically
    a few hundred rows rather than 11M -- so PLINK's parse cost is negligible
    and only the LD-reference load remains.
    """
    sig = df[pd.to_numeric(df["P"], errors="coerce") < p1]
    if sig.empty:
        return []

    if ld_ref and Path(f"{ld_ref}.bed").exists() and shutil.which("plink"):
        assoc = workdir / f"assoc_{tag}.tsv"
        sig[["SNP", "CHR", "POS", "P"]].to_csv(assoc, sep="\t", index=False)
        out = workdir / f"clumped_{tag}"
        # --threads/--memory are not optional here: workers run concurrently and
        # PLINK's defaults are "all cores" and "half of system RAM" each.
        res = run(f"plink --bfile {ld_ref!r} --clump {str(assoc)!r} "
                  f"--clump-snp-field SNP --clump-field P "
                  f"--clump-p1 {p1} --clump-p2 {p1} --clump-r2 {r2} --clump-kb {kb} "
                  f"--threads 1 --memory {plink_mb} "
                  f"--out {str(out)!r}", log)
        clumped = Path(f"{out}.clumped")
        if clumped.exists():
            idx = [ln.split()[2] for ln in clumped.read_text().splitlines()[1:] if ln.strip()]
            return [x for x in idx if x]
        with log.open("a") as fh:
            fh.write(f"plink --clump gave no output for {tag} (rc={res.returncode})\n")
        return []

    return distance_prune(sig, kb)


def distance_prune(sig: pd.DataFrame, kb: int) -> list[str]:
    """Greedy: strongest first, keep if no kept index within +/-kb on the chromosome."""
    d = sig.copy()
    d["P"] = pd.to_numeric(d["P"], errors="coerce")
    d = d.dropna(subset=["P"]).sort_values("P")
    kept: list[tuple[str, int]] = []
    out: list[str] = []
    w = kb * 1000
    for snp, chrom, pos in zip(d["SNP"], d["CHR"].astype(str), d["POS"].astype(int)):
        if any(c == chrom and abs(pos - p) <= w for c, p in kept):
            continue
        kept.append((chrom, pos))
        out.append(snp)
    return out


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
                "var_z_p": np.nan, "ks_p": np.nan, "lambda": np.nan,
                "n_p_lt_05": 0, "frac_p_lt_05": np.nan}
    p = 2.0 * stats.norm.sf(np.abs(z))
    v = float(np.var(z, ddof=1))
    m = z.size
    # Formal one-sided test of var(z) > 1. The null here is ANALYTIC, not
    # estimated: disjoint samples give independent estimates, so var(z)=1
    # exactly. (m-1)*var(z) is therefore chi-square on m-1 df under H0.
    var_p = float(stats.chi2.sf((m - 1) * v, m - 1))
    return {
        "set": label,
        "n": int(m),
        "mean_z": float(np.mean(z)),
        "var_z": v,
        "var_z_p": var_p,
        "ks_p": float(stats.kstest(z, "norm").pvalue),
        "lambda": float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1)),
        "n_p_lt_05": int((p < 0.05).sum()),
        "frac_p_lt_05": float((p < 0.05).mean()),
    }


def variance_ratio_test(rows: list[dict], num: str, den: str) -> dict | None:
    """F test: is one set's between-stratum heterogeneity greater than another's?

    This is the sharpest form of the targeted hypothesis. Asking whether
    smoking-locus pQTLs are modified at all (var > 1) conflates a real effect
    with anything that inflates z generally; asking whether they are modified
    MORE than a random slice of the proteome isolates what is specific to those
    loci, because any shared technical inflation cancels in the ratio.
    """
    a = next((r for r in rows if r["set"] == num), None)
    b = next((r for r in rows if r["set"] == den), None)
    if not a or not b or a["n"] < 3 or b["n"] < 3:
        return None
    if not (np.isfinite(a["var_z"]) and np.isfinite(b["var_z"])) or b["var_z"] <= 0:
        return None
    F = a["var_z"] / b["var_z"]
    return {
        "numerator": num, "denominator": den,
        "n_num": a["n"], "n_den": b["n"],
        "var_num": a["var_z"], "var_den": b["var_z"],
        "F": float(F),
        "p": float(stats.f.sf(F, a["n"] - 1, b["n"] - 1)),
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


def process_protein(protein: str, cfg: dict) -> dict:
    """All work for one protein. Top-level so multiprocessing can pickle it.

    Three streaming passes, one per file, down from six: the reference pass
    yields both the clump input and the reference rows, and each stratum pass
    yields both the index-variant rows and that stratum's own significant
    variants for the secondary tier.
    """
    a, b, ref = cfg["a"], cfg["b"], cfg["ref"]
    pdir = cfg["staged"] / protein
    files = {g: pdir / f"{g}.tsv.gz" for g in (a, b, ref)}

    missing = [g for g, f in files.items() if not f.exists()]
    if missing:
        return {"protein": protein, "loci": None, "skip": f"missing {','.join(missing)}", "msg": ""}
    for g, f in files.items():
        absent = [c for c in REQUIRED if c not in read_header(f)]
        if absent:
            return {"protein": protein, "loci": None,
                    "skip": f"{g} lacks {','.join(absent)}", "msg": ""}

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)

        # Pass 1 -- reference. Everything below the threshold; the index set is
        # necessarily a subset, so no second read is needed.
        dref_all = scan_sumstats(files[ref], work, "ref", cfg["p_discovery"], None, cfg["hla"])
        index = clump_frame(dref_all, work, "ref", cfg["ld_ref"], cfg["p_discovery"],
                            cfg["clump_r2"], cfg["clump_kb"], cfg["log"],
                            cfg["plink_mb"])
        ids = set(index)

        # Passes 2 and 3 -- strata. Each returns rows AT the index variants
        # (needed whatever their p-value there) and rows significant in that
        # stratum (the secondary tier), in a single read.
        dstrata = {g: scan_sumstats(files[g], work, g, cfg["p_discovery"], ids, cfg["hla"])
                   for g in (a, b)}

        # Secondary tier: significant within a stratum but not in the reference.
        # Distance-pruned rather than LD-clumped -- it is exploratory, and a
        # PLINK call per stratum would triple the LD-reference loads.
        sec: list[str] = []
        for g in (a, b):
            d = dstrata[g]
            sig = d[pd.to_numeric(d["P"], errors="coerce") < cfg["p_discovery"]]
            sec += distance_prune(sig, cfg["clump_kb"])
        sec_only = sorted(set(sec) - ids)

        tiers = {s: "primary" for s in index}
        tiers.update({s: "secondary" for s in sec_only})
        if not tiers:
            return {"protein": protein, "loci": None,
                    "skip": "no variant reached the threshold", "msg": ""}

        want = set(tiers)
        # Secondary variants are absent from the reference frame when they sit
        # above its threshold, and from the other stratum's frame; fetch those
        # few by ID. Only happens when a secondary tier exists.
        need_ref = want - set(dref_all["SNP"])
        if need_ref:
            extra = scan_sumstats(files[ref], work, "ref2", 0.0, need_ref, cfg["hla"])
            dref_all = pd.concat([dref_all, extra], ignore_index=True)
        for g in (a, b):
            need = want - set(dstrata[g]["SNP"])
            if need:
                extra = scan_sumstats(files[g], work, f"{g}2", 0.0, need, cfg["hla"])
                dstrata[g] = pd.concat([dstrata[g], extra], ignore_index=True)

        dref = dref_all[dref_all["SNP"].isin(want)]
        if dref.empty:
            return {"protein": protein, "loci": None,
                    "skip": "index variants absent from the reference", "msg": ""}

        da = harmonise(dref, dstrata[a][dstrata[a]["SNP"].isin(want)], f"_{a}")
        db = harmonise(dref, dstrata[b][dstrata[b]["SNP"].isin(want)], f"_{b}")

        base = dref.rename(columns={c: f"{c}_{ref}" for c in
                                    ("BETA", "SE", "P", "N", "EAF", "HET_ISQ", "ID")})
        m = base.merge(da, on="SNP", how="inner").merge(db, on="SNP", how="inner")
        if m.empty:
            return {"protein": protein, "loci": None,
                    "skip": "no variant aligned across all three files", "msg": ""}

        m = compare_pair(m, a, b)
        m = qc_flags(m, a, b, cfg["max_het_isq"], cfg["max_eaf_diff"])
        m = annotate_cis(m, protein, cfg["gene_map"], cfg["cis_window_kb"])
        m.insert(0, "PROTEIN", protein)
        m["TIER"] = m["SNP"].map(tiers)
        m["CONTROL"] = protein in cfg["controls"]

        n_pri = int((m["TIER"] == "primary").sum())
        n_cis = int((m["CIS"] == "cis").sum())
        msg = (f"  {protein:<12} {n_pri:>4} primary, {len(m)-n_pri:>3} secondary"
               f"  ({n_cis} cis, {int((m['CIS']=='trans').sum())} trans)")
        return {"protein": protein, "loci": m, "skip": None, "msg": msg}



# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def default_jobs() -> int:
    """Workers default to the ALLOCATION, not the node.

    os.cpu_count() reports the machine's cores even inside a 2-CPU cgroup, so
    under sbatch the old default (cpu_count-1) asked for tens of workers that
    Slurm then squeezed onto the allocated CPUs.
    """
    n = os.environ.get("SLURM_CPUS_PER_TASK")
    if n and n.isdigit() and int(n) > 0:
        return int(n)
    return max(1, (os.cpu_count() or 2) - 1)


def plink_memory_mb(requested: int, jobs: int) -> int:
    """MB to hand PLINK per worker.

    PLINK1.9 with no --memory claims half of TOTAL system RAM. Sixteen of those
    inside one 16G cgroup is an OOM, so the budget is always stated explicitly.
    """
    if requested > 0:
        return requested
    for var in ("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"):
        v = os.environ.get(var, "")
        if v.isdigit() and int(v) > 0:
            total = int(v) * (default_jobs() if var == "SLURM_MEM_PER_CPU" else 1)
            return max(256, int(total * 0.8) // max(1, jobs))
    return 2000


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
    ap.add_argument("--jobs", "-j", type=int, default=default_jobs(),
                    help="parallel workers; proteins are independent. Defaults to "
                         "SLURM_CPUS_PER_TASK under sbatch, else cpu_count-1.")
    ap.add_argument("--plink-memory", type=int, default=0,
                    help="MB per PLINK call. 0 = derive from the SLURM memory "
                         "allocation divided by --jobs, else 2000.")
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
    plink_mb = plink_memory_mb(args.plink_memory, args.jobs)
    print(f"{len(proteins)} protein(s); {len(controls)} flagged as controls")
    print(f"workers: {args.jobs}   PLINK: 1 thread, {plink_mb} MB each\n")

    cfg = dict(a=a, b=b, ref=ref, staged=staged, log=log, gene_map=gene_map,
               controls=controls,
               p_discovery=args.p_discovery, ld_ref=args.ld_ref,
               clump_r2=args.clump_r2, clump_kb=args.clump_kb, hla=args.hla,
               cis_window_kb=args.cis_window_kb,
               max_het_isq=args.max_het_isq, max_eaf_diff=args.max_eaf_diff,
               plink_mb=plink_mb)

    if args.jobs == 1:
        results = [process_protein(prot, cfg) for prot in proteins]
    else:
        # Proteins are independent, and the work is dominated by decompressing
        # one file per stratum, so this scales close to linearly until the
        # filesystem saturates.
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(args.jobs) as pool:
            results = pool.starmap(process_protein, [(prot, cfg) for prot in proteins],
                                   chunksize=1)

    all_loci = [r["loci"] for r in results if r["loci"] is not None]
    skipped = [(r["protein"], r["skip"]) for r in results if r["skip"]]
    for r in results:
        if r["msg"]:
            print(r["msg"])

    if not all_loci:
        for p, why in skipped:
            print(f"  {p:<12} skipped: {why}")
        sys.exit("\nNo protein yielded index variants; nothing to compare.")

    loci = pd.concat(all_loci, ignore_index=True)

    # Multiple testing over PRIMARY, QC-passing, non-control variants only:
    # those are the actual hypothesis tests.
    testable = loci["QC_PASS"] & (loci["TIER"] == "primary") & (~loci["CONTROL"])
    n_tests = int(testable.sum())
    # Correct within cis and within trans, matching how they are reported. A
    # single pooled denominator would penalise a cis finding for the number of
    # trans tests, which answer a different question.
    loci["P_DIFF_BONF"] = np.nan
    for grp in ("cis", "trans", "NA"):
        sel = testable & (loci["CIS"] == grp)
        k = int(sel.sum())
        if k:
            loci.loc[sel, "P_DIFF_BONF"] = np.minimum(loci.loc[sel, "P_DIFF"] * k, 1.0)
    loci["N_TESTS_IN_GROUP"] = np.where(
        testable, loci["CIS"].map(
            {g: int((testable & (loci["CIS"] == g)).sum()) for g in loci["CIS"].unique()}
        ), np.nan)
    loci["SIGNIF_DIFF"] = loci["P_DIFF_BONF"] < 0.05

    loci = loci.sort_values(["PROTEIN", "P_DIFF"])

    front = ["PROTEIN", "SNP", "CHR", "POS", "EA", "OA", "GENE", "CIS",
             "TIER", "CONTROL"]
    stats_cols = ["DELTA", "SE_DELTA", "Z_DIFF", "P_DIFF", "P_DIFF_BONF",
                  "N_TESTS_IN_GROUP", "SIGNIF_DIFF", "SAME_DIRECTION", "Q", "I2"]
    per_group = [f"{m}_{g}" for g in (a, b, ref)
                 for m in ("BETA", "SE", "P", "N", "EAF", "HET_ISQ")]
    qc = ["EAF_DIFF", "MAX_HET_ISQ", "QC_EAF", "QC_HET", "QC_PASS"]
    order = [c for c in front + stats_cols + per_group + qc if c in loci.columns]
    order += [c for c in loci.columns if c not in order]
    loci = loci[order]

    loci.to_csv(outdir / "pqtl_comparison.tsv", sep="\t", index=False, na_rep="NA")

    # Global tests: the primary question, and the control set as empirical null.
    # cis and trans are reported separately, not pooled. They are different
    # biology -- cis acts near the gene and is largely context-independent,
    # whereas trans acts through regulatory intermediates and is where an
    # exposure could plausibly intervene -- and they have very different power,
    # since cis effects are several times larger. Pooling averages any
    # trans-specific signal away. Calibration is unaffected: var(z)=1 under the
    # null regardless of effect size, so each stratum of the table has the same
    # reference point.
    cis_m   = loci["CIS"] == "cis"
    trans_m = loci["CIS"] == "trans"
    ctrl    = loci["QC_PASS"] & loci["CONTROL"]
    secnd   = loci["QC_PASS"] & (loci["TIER"] == "secondary")

    rows = [
        global_test(loci.loc[testable, "Z_DIFF"].values, "primary_targeted"),
        global_test(loci.loc[testable & cis_m, "Z_DIFF"].values, "primary_targeted_cis"),
        global_test(loci.loc[testable & trans_m, "Z_DIFF"].values, "primary_targeted_trans"),
        global_test(loci.loc[ctrl, "Z_DIFF"].values, "primary_random_proteome"),
        global_test(loci.loc[ctrl & cis_m, "Z_DIFF"].values, "primary_random_proteome_cis"),
        global_test(loci.loc[ctrl & trans_m, "Z_DIFF"].values, "primary_random_proteome_trans"),
        global_test(loci.loc[secnd, "Z_DIFF"].values, "secondary_selection_biased"),
    ]
    gt = pd.DataFrame(rows)
    gt.to_csv(outdir / "global_test.tsv", sep="\t", index=False, na_rep="NA")

    # Is the targeted set MORE heterogeneous than a random slice of the
    # proteome? Any inflation common to both cancels in the ratio, so this
    # isolates what is specific to the smoking loci.
    ratios = [r for r in (
        variance_ratio_test(rows, "primary_targeted", "primary_random_proteome"),
        variance_ratio_test(rows, "primary_targeted_cis", "primary_random_proteome_cis"),
        variance_ratio_test(rows, "primary_targeted_trans", "primary_random_proteome_trans"),
    ) if r]
    if ratios:
        pd.DataFrame(ratios).to_csv(outdir / "variance_ratio.tsv", sep="\t",
                                    index=False, na_rep="NA")

    plotted = make_plots(loci, outdir, a, b)

    # ---- report ----------------------------------------------------------
    print(f"\n{'='*78}\nGLOBAL TESTS -- var(z) = 1 under the null (analytic, not estimated)")
    print(f"{'set':<34}{'n':>5}{'mean z':>9}{'var z':>8}{'P(var>1)':>10}{'p<.05':>8}")
    for r in rows:
        if r["n"] == 0:
            continue
        print(f"{r['set']:<34}{r['n']:>5}{r['mean_z']:>9.3f}{r['var_z']:>8.2f}"
              f"{r['var_z_p']:>10.3g}{r['frac_p_lt_05']:>8.2f}")

    print(f"\n{'-'*78}\nTHE THREE QUESTIONS")
    def get(lbl):
        return next((r for r in rows if r["set"] == lbl and r["n"] >= 3), None)

    def verdict(pv: float) -> str:
        """Graded rather than binary. var(z) is a single test on a quantity with
        a wide sampling distribution -- at n=290 a true null produces var(z)=1.16
        about 3.5% of the time -- so a p just under 0.05 is weak evidence and
        should not be reported as a finding."""
        if not np.isfinite(pv):
            return "not computable"
        if pv < 0.001:
            return "YES"
        if pv < 0.05:
            return "suggestive only"
        return "no evidence"

    rnd, tgt = get("primary_random_proteome"), get("primary_targeted")
    if rnd:
        print(f"  1. Do pQTLs differ by stratum GENERALLY?        {verdict(rnd['var_z_p'])}")
        print(f"     random proteome: n={rnd['n']}, var(z)={rnd['var_z']:.2f}, P={rnd['var_z_p']:.3g}")
    else:
        print("  1. Do pQTLs differ GENERALLY?                   no random set "
              "(pass --annot)")
    if tgt:
        print(f"  2. Do pQTLs at the TARGETED loci differ?        {verdict(tgt['var_z_p'])}")
        print(f"     targeted: n={tgt['n']}, var(z)={tgt['var_z']:.2f}, P={tgt['var_z_p']:.3g}")
    if ratios:
        r0 = ratios[0]
        print(f"  3. Do the TARGETED loci differ MORE than average?  {verdict(r0['p'])}")
        print(f"     F = {r0['var_num']:.2f}/{r0['var_den']:.2f} = {r0['F']:.2f}"
              f" on ({r0['n_num']-1},{r0['n_den']-1}) df, P={r0['p']:.3g}")
    else:
        print("  3. Targeted vs random ratio: not computable (need both sets)")

    print("\n  Question 1 uses the random proteome as a DISCOVERY set, so it cannot")
    print("  simultaneously serve as an empirical null. It is interpretable only")
    print("  because the null is analytic and the QC gates below passed; a")
    print("  permutation of stratum labels would be needed for an empirical one.")

    n_fail = int((~loci["QC_PASS"]).sum())
    print(f"\n{'='*70}\nPER-LOCUS")
    print(f"  index variants            : {len(loci)}")
    print(f"  primary / secondary       : {int((loci['TIER']=='primary').sum())}"
          f" / {int((loci['TIER']=='secondary').sum())}")
    print(f"  failed QC (EAF or I2)     : {n_fail}")
    print(f"  hypothesis tests          : {n_tests}"
          f"  (cis {int((testable & cis_m).sum())}, trans {int((testable & trans_m).sum())})")
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
