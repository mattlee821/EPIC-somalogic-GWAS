#!/usr/bin/env python3
"""Compare the observed statistics against their permutation distributions.

Two observed values are reported for each statistic and they answer different
questions:

  OLS      the same engine that produced the permutations, on the TRUE labels.
           This is the like-for-like comparison and the empirical P comes from it.
  PIPELINE what compare_strata.py reported, from REGENIE + METAL. Shown so any
           gap between the two estimators is visible rather than assumed away.

A large OLS-vs-PIPELINE gap does not invalidate the permutation; it means the
OLS approximation (no LOCO term, one-pass IVW) shifts that statistic, and the
empirical P should be read against the OLS column.
"""
from __future__ import annotations
import argparse, csv, glob, math
from pathlib import Path
import numpy as np


def load(pat):
    rows = []
    for f in sorted(glob.glob(pat)):
        with open(f) as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    return rows


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm-dir", required=True)
    ap.add_argument("--comparison-dir", default="")
    a = ap.parse_args()
    pdir = Path(a.perm_dir)

    perm = load(str(pdir / "perm_*.tsv"))
    if not perm:
        raise SystemExit(f"no perm_*.tsv under {pdir}")
    obs_rows = load(str(pdir / "observed_ols.tsv"))
    obs = obs_rows[0] if obs_rows else None

    # what the pipeline reported, for the side-by-side
    pipe = {}
    cdir = Path(a.comparison_dir) if a.comparison_dir else pdir.parent / "comparison"
    gt = cdir / "global_test.tsv"
    if gt.exists():
        for r in load(str(gt)):
            s = r["set"]
            if s == "primary_random_proteome_dedup": pipe["var_z_random_dedup"] = fnum(r["var_z"])
            if s == "primary_targeted_dedup":        pipe["var_z_targeted_dedup"] = fnum(r["var_z"])
            if s == "primary_random_proteome":       pipe["var_z_random"] = fnum(r["var_z"])
            if s == "primary_targeted":              pipe["var_z_targeted"] = fnum(r["var_z"])
    om = cdir / "oriented_mean.tsv"
    if om.exists():
        for r in load(str(om)):
            if r["set"] == "primary_random_proteome":
                pipe["mean_oriented_random"] = fnum(r["mean_z_oriented"])
            if r["set"] == "primary_targeted":
                pipe["mean_oriented_targeted"] = fnum(r["mean_z_oriented"])

    stats = ["var_z_random_dedup", "var_z_targeted_dedup",
             "mean_oriented_random", "mean_oriented_targeted", "F_ratio"]
    stats += [c for c in perm[0] if c.startswith("var_z_maf_")]

    n = len(perm)
    print(f"permutations: {n}")
    if n < 200:
        print("  WARNING: fewer than 200 permutations -- the 2.5/97.5 percentiles")
        print("           are estimated from a handful of values.")
    print()
    print(f"{'statistic':<26}{'OLS obs':>9}{'pipeline':>10}{'null mean':>11}"
          f"{'null 2.5%':>11}{'null 97.5%':>11}{'emp. P':>9}")
    print("-" * 87)
    out = []
    for s in stats:
        v = np.array([fnum(r.get(s)) for r in perm])
        v = v[np.isfinite(v)]
        if v.size < 20:
            continue
        lo, hi = np.percentile(v, [2.5, 97.5])
        o = fnum(obs.get(s)) if obs else np.nan
        # one-sided for variances and the F ratio (inflation is the alternative),
        # two-sided for the oriented means (direction is not specified a priori)
        if np.isfinite(o):
            if s.startswith("mean_oriented"):
                p = (1 + np.sum(np.abs(v - v.mean()) >= abs(o - v.mean()))) / (1 + v.size)
            else:
                p = (1 + np.sum(v >= o)) / (1 + v.size)
        else:
            p = np.nan
        pp = pipe.get(s, np.nan)
        f3 = lambda x: f"{x:>9.3f}" if np.isfinite(x) else f"{'--':>9}"
        g3 = lambda x: f"{x:>10.3f}" if np.isfinite(x) else f"{'--':>10}"
        ps = f"{p:>9.4f}" if np.isfinite(p) else f"{'--':>9}"
        print(f"{s:<26}{f3(o)}{g3(pp)}{v.mean():>11.3f}{lo:>11.3f}{hi:>11.3f}{ps}")
        out.append(dict(statistic=s, observed_ols=o, observed_pipeline=pp,
                        null_mean=v.mean(), null_sd=v.std(ddof=1),
                        null_lo=lo, null_hi=hi, emp_p=p, n_perm=int(v.size)))

    if not out:
        raise SystemExit("no statistic had enough finite permutations to summarise")
    with (pdir / "permutation_null.tsv").open("w") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0]))
        w.writeheader()
        for r in out:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})

    print("\n" + "=" * 87)
    ref = next((r for r in out if r["statistic"] == "var_z_random_dedup"), None)
    if ref:
        print(f"CALIBRATION  the random proteome under permutation centres on "
              f"{ref['null_mean']:.3f}")
        print("             (the analytic null asserts 1.000)")
    print("\nNOTE  the refit is OLS without the REGENIE LOCO term and uses one-pass")
    print("      IVW. Read every empirical P against the OLS column, not the")
    print("      pipeline column.")
    print(f"\nwrote {pdir/'permutation_null.tsv'}")


if __name__ == "__main__":
    main()
