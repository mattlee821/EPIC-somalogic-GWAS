#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from config import ensure_string_ids, make_safe_names, read_table_auto


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cov_file", required=True)
    parser.add_argument("--pcs_file", required=True)
    parser.add_argument("--eval_file", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def load_pcs_table(pcs_file: str) -> pd.DataFrame:
    pcs = read_table_auto(pcs_file)
    if "#FID" in pcs.columns and "FID" not in pcs.columns:
        pcs = pcs.rename(columns={"#FID": "FID"})

    if {"FID", "IID"}.issubset(pcs.columns):
        return pcs

    pcs = pd.read_csv(pcs_file, sep=r"\s+", engine="python", header=None)
    if pcs.shape[1] < 3:
        raise ValueError("PC file has too few columns (expected FID IID PC1 ...).")

    pcs.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pcs.shape[1] - 1)]
    return pcs


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    covariates = read_table_auto(args.cov_file)
    pcs = load_pcs_table(args.pcs_file)

    covariates = ensure_string_ids(covariates)
    pcs = ensure_string_ids(pcs)

    merged = covariates.merge(pcs, on=["FID", "IID"], how="left", sort=False)

    other_columns = [col for col in merged.columns if col not in {"FID", "IID"}]
    if other_columns:
        merged = merged.rename(columns=dict(zip(other_columns, make_safe_names(other_columns))))

    merged.to_csv(outdir / "covariates.cov", sep="\t", index=False, na_rep="NA")

if __name__ == "__main__":
    main()
