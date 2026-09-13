#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from config import read_table_auto


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics_files", default="")
    parser.add_argument("--metrics_list", default="")
    # Consolidated per-protein LDSC / heritability rows (each a header + N rows,
    # already merged by collectFile). Empty/placeholder => that side is skipped.
    parser.add_argument("--ldsc_file", default="")
    parser.add_argument("--herit_file", default="")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def _read_metrics_table(path: str) -> pd.DataFrame | None:
    """Read a metrics table if it exists and has at least one data row."""
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    df = read_table_auto(str(p))
    if df is None or df.empty or "protein" not in df.columns:
        return None
    return df


def write_ldsc_herit_summary(ldsc_file: str, herit_file: str, outdir: Path) -> None:
    """Join the LDSC and heritability per-protein metrics into one table.

    Both tables share intercept/lambda column names but mean different things
    (LDSC = full-genome; heritability = polygenic remainder), so every non-key
    column is namespaced with an ldsc_/herit_ prefix before the outer join.

    The join key is (protein, group). Keyed on protein alone, a run with two
    strata cross-joins LDSC's two rows against heritability's two rows and
    yields four meaningless rows per protein, with no way to tell which stratum
    any figure belongs to.
    """
    ldsc = _read_metrics_table(ldsc_file)
    herit = _read_metrics_table(herit_file)
    if ldsc is None and herit is None:
        return

    # Join on whichever key columns both sides actually carry, so rows written
    # before 'group' existed still merge instead of raising.
    keys = [k for k in ("protein", "group")
            if (ldsc is None or k in ldsc.columns) and (herit is None or k in herit.columns)]

    def _prefix(df: pd.DataFrame, pre: str) -> pd.DataFrame:
        return df.rename(columns={c: f"{pre}{c}" for c in df.columns if c not in keys})

    if ldsc is not None and herit is not None:
        merged = pd.merge(_prefix(ldsc, "ldsc_"), _prefix(herit, "herit_"),
                          on=keys, how="outer")
    elif ldsc is not None:
        merged = _prefix(ldsc, "ldsc_")
    else:
        merged = _prefix(herit, "herit_")

    merged = merged.sort_values(keys)
    merged.to_csv(outdir / "metrics_ldsc-herit.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Study-specific GWAS-QC metrics -> all_metrics.tsv
    metrics_paths: list[str] = []
    if args.metrics_list and Path(args.metrics_list).exists():
        metrics_paths = [line.strip() for line in Path(args.metrics_list).read_text().splitlines() if line.strip()]
    elif args.metrics_files:
        metrics_paths = [part.strip() for part in args.metrics_files.split(",") if part.strip()]

    if metrics_paths:
        metrics_frames = [read_table_auto(path) for path in metrics_paths]
        metrics_all = pd.concat(metrics_frames, ignore_index=True, sort=False)
        drop_columns = ["flag_raw", "flag_filtered", "info_score"]
        present = [column for column in drop_columns if column in metrics_all.columns]
        if present:
            metrics_all = metrics_all.drop(columns=present)
        metrics_all.to_csv(outdir / "all_metrics.tsv", sep="\t", index=False)

    # 2. LDSC + heritability metrics -> metrics_ldsc-herit.tsv
    write_ldsc_herit_summary(args.ldsc_file, args.herit_file, outdir)


if __name__ == "__main__":
    main()
