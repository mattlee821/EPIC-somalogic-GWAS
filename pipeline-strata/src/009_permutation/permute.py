#!/usr/bin/env python3
"""Permutation null for the stratified pQTL comparison (METHODS.md 8.3).

Permutes smoking labels WITHIN study x case-status cell, refits at the fixed
index variants, and recomputes every statistic the comparison reports. The index
set is held at the observed one, so this is conditional inference: the null is
"no stratum difference GIVEN this index set".

Two approximations, both declared rather than hidden:

  * The refit is ordinary least squares. The real analysis carries a REGENIE
    LOCO polygenic term, which cannot be refitted per permutation (step 1 would
    have to run 16 times per permutation). The LOCO term is a per-individual
    offset independent of the tested genotype, so it largely cancels in a
    BETWEEN-STRATUM difference -- but "largely" is not "exactly".
  * Meta-analysis is inverse-variance weighted over studies, matching METAL's
    SCHEME STDERR, but computed here in one pass rather than the pipeline's
    two-stage within-study-then-across-study route. Fixed-effect IVW is
    associative, so these agree.

Speed comes from Frisch-Waugh: residualise Y and G on the covariates once per
study x stratum, after which every beta is a single dot product instead of a
fitted regression. That turns ~54,000 regressions per permutation into 16
residualisation passes.
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- loading
def read_table(path: Path, sep="\t"):
    """Header + float/str columns, without a pandas dependency."""
    with path.open() as fh:
        head = fh.readline().rstrip("\n").split(sep)
        rows = [ln.rstrip("\n").split(sep) for ln in fh if ln.strip()]
    return head, rows


def load_study(study: str, pdir: Path, want_pairs, want_prot):
    """Covariates, phenotypes and harmonised dosages for one study."""
    g = pdir / "geno"
    raw_p, pvar_p = g / f"{study}.raw", g / f"{study}.pvar"
    cov_p = pdir / "cov" / f"{study}.cov"
    phe_p = pdir / "pheno" / f"{study}.pheno"
    lab_p = pdir / "labels" / f"{study}.tsv"
    for p in (raw_p, pvar_p, cov_p, phe_p, lab_p):
        if not p.exists():
            raise FileNotFoundError(p)

    # --- variants: position -> the column whose alleles match the sumstats ---
    # A .pvar opens with '##' metadata lines; the real header is the '#CHROM'
    # line. Taking line 1 as the header silently yields a metadata line.
    pv_head, pv_rows = None, []
    with pvar_p.open() as fh:
        for ln in fh:
            if ln.startswith("##"):
                continue
            f = ln.rstrip("\n").split("\t")
            if pv_head is None and (ln.startswith("#") or f[0].upper() in ("CHROM", "#CHROM")):
                pv_head = [h.lstrip("#").upper() for h in f]
                continue
            if pv_head is None:                      # headerless: assume .bim order
                pv_head = ["CHROM", "ID", "CM", "POS", "ALT", "REF"][:len(f)]
            pv_rows.append(f)
    need = ("CHROM", "POS", "ID", "REF", "ALT")
    missing = [c for c in need if c not in pv_head]
    if missing:
        raise RuntimeError(f"{pvar_p}: header lacks {missing}; found {pv_head}")
    ix = {c: pv_head.index(c) for c in need}
    by_pos: dict[str, list] = {}
    for r in pv_rows:
        key = f"{r[ix['CHROM']].replace('chr','')}:{r[ix['POS']]}"
        by_pos.setdefault(key, []).append(
            (r[ix["ID"]], r[ix["REF"]].upper(), r[ix["ALT"]].upper()))

    raw_head, raw_rows = read_table(raw_p)
    # plink2 --export A writes "<ID>_<countedallele>"; counted allele is ALT.
    col_of = {}
    for j, c in enumerate(raw_head[6:], start=6):
        vid = c.rsplit("_", 1)[0]
        col_of.setdefault(vid, j)

    snp_col, snp_sign, dropped = {}, {}, 0
    for snp, ea, oa in want_pairs:
        if snp in snp_col:
            continue
        hit = None
        for vid, ref, alt in by_pos.get(snp, []):
            if {ref, alt} == {ea.upper(), oa.upper()} and vid in col_of:
                hit = (col_of[vid], 1.0 if alt == ea.upper() else -1.0)
                break
        if hit is None:
            dropped += 1
            continue
        snp_col[snp], snp_sign[snp] = hit

    iid = [r[1] for r in raw_rows]
    G = np.empty((len(raw_rows), len(snp_col)), dtype=np.float64)
    order = list(snp_col)
    for k, snp in enumerate(order):
        j, s = snp_col[snp], snp_sign[snp]
        colv = np.array([(np.nan if r[j] in ("NA", "") else float(r[j])) for r in raw_rows])
        # flip so the dosage always counts the summary-stat effect allele
        G[:, k] = colv if s > 0 else 2.0 - colv
    gmean = np.nanmean(G, axis=0)
    miss = np.isnan(G)
    if miss.any():
        G[miss] = np.take(gmean, np.where(miss)[1])
    gidx = {snp: k for k, snp in enumerate(order)}

    # --- covariates ---
    ch, cr = read_table(cov_p)
    cmap = {r[1]: r for r in cr}
    cov_cols = [i for i, c in enumerate(ch) if c not in ("FID", "IID")]
    # --- phenotypes ---
    ph, pr = read_table(phe_p)
    pmap = {r[1]: r for r in pr}
    pcol = {c: i for i, c in enumerate(ph)}
    prot = [p for p in want_prot if p in pcol]
    pidx = {p: k for k, p in enumerate(prot)}
    # --- labels ---
    lh, lr = read_table(lab_p)
    lab = {r[0]: (r[1], r[2]) for r in lr}          # iid -> (SMOKE, PHENO cell)

    keep = [i for i, s in enumerate(iid) if s in cmap and s in pmap and s in lab]
    if not keep:
        raise RuntimeError(f"{study}: no samples shared between geno/cov/pheno/labels")
    G = G[keep]
    ids = [iid[i] for i in keep]
    X = np.array([[float(cmap[s][i]) for i in cov_cols] for s in ids])
    X = np.column_stack([np.ones(len(ids)), X])
    Y = np.array([[float(pmap[s][pcol[p]]) for p in prot] for s in ids])
    smoke = np.array([lab[s][0] for s in ids])
    cell = np.array([lab[s][1] for s in ids])
    return dict(study=study, ids=ids, X=X, Y=Y, G=G, gidx=gidx, pidx=pidx,
                smoke=smoke, cell=cell, dropped=dropped, gmean=gmean, order=order)


# ---------------------------------------------------------------- the fit
def stratum_betas(X, Y, G, rows, pi, gi):
    """OLS beta/se for each (protein, variant) pair, via Frisch-Waugh."""
    Xs, Ys, Gs = X[rows], Y[rows], G[rows]
    n, p = Xs.shape
    dof = n - p - 1
    if dof < 20:
        return None, None
    Q, _ = np.linalg.qr(Xs)
    Yr = Ys - Q @ (Q.T @ Ys)
    Gr = Gs - Q @ (Q.T @ Gs)
    Gp, Yp = Gr[:, gi], Yr[:, pi]
    gg = np.einsum("ij,ij->j", Gp, Gp)
    gg[gg <= 1e-8] = np.nan
    beta = np.einsum("ij,ij->j", Gp, Yp) / gg
    yy = np.einsum("ij,ij->j", Yp, Yp)
    rss = np.maximum(yy - beta ** 2 * gg, 1e-12)
    se = np.sqrt(rss / dof / gg)
    return beta, se


def ivw(betas, ses):
    """Inverse-variance weighted fixed-effect meta over studies."""
    w = np.zeros_like(betas[0]); bw = np.zeros_like(betas[0])
    for b, s in zip(betas, ses):
        ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
        wi = np.where(ok, 1.0 / np.where(ok, s, 1.0) ** 2, 0.0)
        bw += np.where(ok, b, 0.0) * wi
        w += wi
    w[w <= 0] = np.nan
    return bw / w, np.sqrt(1.0 / w)


# ---------------------------------------------------------------- stats
def summarise(z, bmeta, targeted, dedup_keep, maf, bands):
    out = {}
    fin = np.isfinite(z)
    def vz(m):
        v = z[m & fin]
        return (float(np.var(v, ddof=1)), int(v.size)) if v.size >= 3 else (np.nan, int(v.size))
    def om(m):
        v = z[m & fin] * np.sign(bmeta[m & fin])
        return float(np.mean(v)) if v.size >= 3 else np.nan
    rnd = ~targeted
    out["var_z_random"], out["n_random"] = vz(rnd)
    out["var_z_targeted"], out["n_targeted"] = vz(targeted)
    out["var_z_random_dedup"], _ = vz(rnd & dedup_keep)
    out["var_z_targeted_dedup"], _ = vz(targeted & dedup_keep)
    out["mean_oriented_random"] = om(rnd)
    out["mean_oriented_targeted"] = om(targeted)
    vt, vr = out["var_z_targeted_dedup"], out["var_z_random_dedup"]
    out["F_ratio"] = vt / vr if (np.isfinite(vt) and np.isfinite(vr) and vr > 0) else np.nan
    for lo, hi in bands:
        m = rnd & (maf > lo) & (maf <= hi)
        out[f"var_z_maf_{lo}_{hi}"], _ = vz(m)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm-dir", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--observed", action="store_true",
                    help="run once on the TRUE labels, no shuffle. Needed so the "
                         "observed statistic is computed by the same OLS engine "
                         "as the permutations -- comparing a REGENIE/METAL "
                         "statistic against an OLS null is not like for like.")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    pdir = Path(a.perm_dir)
    ph, pr = read_table(pdir / "pairs.tsv")
    pairs = [(r[0], r[1], r[2], r[3]) for r in pr]
    want_pairs = [(s, ea, oa) for _, s, ea, oa in pairs]
    want_prot = sorted({p for p, _, _, _ in pairs})

    ah, ar = read_table(Path(a.annot))
    aix = {c: ah.index(c) for c in ("seqid", "source")}
    src = {r[aix["seqid"]]: r[aix["source"]] for r in ar}
    targeted = np.array([src.get(p, "") == "targeted" for p, _, _, _ in pairs])

    studies = sorted(p.stem for p in (pdir / "geno").glob("*.raw"))
    log(f"studies: {len(studies)}   pairs: {len(pairs)}   proteins: {len(want_prot)}")
    S = []
    for s in studies:
        d = load_study(s, pdir, want_pairs, want_prot)
        log(f"  {s:<10} n={len(d['ids']):>5}  variants={d['G'].shape[1]:>5}  "
            f"proteins={d['Y'].shape[1]:>4}  unmatched variants={d['dropped']}")
        S.append(d)

    # per-study index vectors for the pair list
    for d in S:
        pi, gi, ok = [], [], []
        for prot, snp, _, _ in pairs:
            if prot in d["pidx"] and snp in d["gidx"]:
                pi.append(d["pidx"][prot]); gi.append(d["gidx"][snp]); ok.append(True)
            else:
                pi.append(0); gi.append(0); ok.append(False)
        d["pi"] = np.array(pi); d["gi"] = np.array(gi); d["ok"] = np.array(ok)

    # allele frequency for the MAF bands, pooled over studies
    num = np.zeros(len(pairs)); den = np.zeros(len(pairs))
    for d in S:
        f = d["G"].mean(axis=0) / 2.0
        num += np.where(d["ok"], f[d["gi"]] * len(d["ids"]), 0.0)
        den += np.where(d["ok"], len(d["ids"]), 0.0)
    eaf = np.divide(num, den, out=np.full(len(pairs), np.nan), where=den > 0)
    maf = np.minimum(eaf, 1 - eaf)
    bands = [(0.01, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5)]

    # one fixed deduplication, applied identically to every permutation
    rng0 = np.random.default_rng(0)
    seen, keep = set(), np.zeros(len(pairs), dtype=bool)
    for i in rng0.permutation(len(pairs)):
        snp = pairs[i][1]
        if snp not in seen:
            seen.add(snp); keep[i] = True
    log(f"dedup keeps {keep.sum()} of {len(pairs)} pairs\n")

    rows = []
    t0 = time.time()
    if a.observed:
        log("OBSERVED run: true labels, no shuffle\n")
    for k in range(a.start, a.start + (1 if a.observed else a.n)):
        rng = np.random.default_rng(a.seed + k)
        be, se_, bn, sn = [], [], [], []
        for d in S:
            sm = d["smoke"].copy()
            if not a.observed:
                # PERMUTE WITHIN study x case-status cell: a global shuffle would
                # erase the study-level stratum imbalance the design enforces and
                # the smoking-disease association, giving a null too narrow.
                for c in np.unique(d["cell"]):
                    m = np.where(d["cell"] == c)[0]
                    sm[m] = rng.permutation(sm[m])
            for arm, B, E in (("EVER", be, se_), ("NEVER", bn, sn)):
                r = np.where(sm == arm)[0]
                b, s = stratum_betas(d["X"], d["Y"], d["G"], r, d["pi"], d["gi"])
                if b is None:
                    continue
                b = np.where(d["ok"], b, np.nan); s = np.where(d["ok"], s, np.nan)
                B.append(b); E.append(s)
        if not be or not bn:
            continue
        b_e, s_e = ivw(be, se_)
        b_n, s_n = ivw(bn, sn)
        sd = np.sqrt(s_e ** 2 + s_n ** 2)
        z = (b_e - b_n) / sd
        wm = 1 / s_e ** 2 + 1 / s_n ** 2
        b_meta = (b_e / s_e ** 2 + b_n / s_n ** 2) / wm
        r = summarise(z, b_meta, targeted, keep, maf, bands)
        r["perm"] = k
        rows.append(r)
        if (k - a.start + 1) % 10 == 0:
            log(f"  {k - a.start + 1}/{a.n} perms  ({time.time()-t0:.0f}s)")

    out = (Path(a.out) if a.out else
           pdir / ("observed_ols.tsv" if a.observed else f"perm_{a.start:05d}.tsv"))
    cols = ["perm"] + [c for c in rows[0] if c != "perm"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")
    log(f"\nwrote {out}  ({len(rows)} permutations, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
