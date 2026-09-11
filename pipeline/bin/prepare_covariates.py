#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from config import ensure_string_ids, read_path_list, read_sample_file, read_table_auto, resolve_column_name, standardize_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--covariate_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--study_id", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--group_column", default="")
    parser.add_argument("--cases_value", default="")
    parser.add_argument("--include_covariates", default="")
    parser.add_argument("--include_proteins", default="")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def binary_to_zero_one(series: pd.Series) -> pd.Series:
    non_missing = series.dropna()
    sorted_values = sorted(non_missing.unique())
    lower_value = sorted_values[0]
    return pd.Series(
        np.where(series.isna(), np.nan, np.where(series == lower_value, 0, 1)),
        index=series.index,
        dtype="float",
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    covariates = standardize_ids(read_table_auto(args.covariate_file))
    samples = read_sample_file(args.sample_file)

    if "IID" not in covariates.columns or "IID" not in samples.columns:
        raise ValueError("Both covariate and sample inputs must contain an IID column.")

    covariates = ensure_string_ids(covariates, ("IID",))
    samples = ensure_string_ids(samples, ("IID",))

    # Restrict to samples assigned to THIS study by the cross-study dedup (003),
    # mirroring prepare_phenotypes. Without it, a person genotyped in several
    # studies (shared Idepic_Bio IID) is matched into every study, undoing dedup.
    if "STUDY" in covariates.columns:
        covariates = covariates[covariates["STUDY"].astype(str) == str(args.study_id)].copy()
        if covariates.empty:
            raise ValueError(f"No covariate rows assigned to study '{args.study_id}'.")

    common_samples = pd.Index(covariates["IID"].dropna().unique()).intersection(samples["IID"].dropna().unique())
    if common_samples.empty:
        raise ValueError("Zero overlap found.")

    covariates = covariates[covariates["IID"].isin(common_samples)].copy()
    # Match the .psam convention (FID == IID) so plink2 --keep / REGENIE join on
    # FID+IID succeed; a hardcoded FID of "0" would match zero samples.
    covariates.insert(0, "FID", covariates["IID"])
    covariates = covariates[["FID", "IID"] + [col for col in covariates.columns if col not in {"FID", "IID"}]]

    if args.include_covariates.strip():
        requested_covariates = [part.strip() for part in args.include_covariates.split(",") if part.strip()]
        if args.group == "combined" and args.group_column and args.group_column not in requested_covariates:
            requested_covariates.append(args.group_column)
        columns_to_keep = ["FID", "IID"] + requested_covariates
        covariates = covariates[[col for col in columns_to_keep if col in covariates.columns]].copy()
    else:
        exclude_columns = {"STUDY", "FID", "IID", "FID_1", "IID_1", "PlateId"}
        proteins_to_exclude = read_path_list(args.include_proteins)
        if proteins_to_exclude is not None:
            exclude_columns.update(proteins_to_exclude)

        numeric_columns = [col for col in covariates.columns if pd.api.types.is_numeric_dtype(covariates[col])]
        keep_columns = [col for col in numeric_columns if col not in exclude_columns]
        covariates = covariates[["FID", "IID"] + keep_columns].copy()

    if args.group and args.group != "combined":
        sample_group_column = resolve_column_name(samples.columns, args.group_column)
        covariate_group_column = resolve_column_name(covariates.columns, args.group_column)

        if sample_group_column:
            group_frame = samples[["IID", sample_group_column]].copy()
            if sample_group_column != args.group_column:
                group_frame = group_frame.rename(columns={sample_group_column: args.group_column})
            covariates = covariates.merge(group_frame, on="IID", how="inner", sort=False)
        elif covariate_group_column:
            if covariate_group_column != args.group_column:
                covariates = covariates.rename(columns={covariate_group_column: args.group_column})
            pass
        else:
            raise ValueError(f"Group column not found in sample or covariate file: {args.group_column}")

        group_values = covariates[args.group_column].astype("string")
        if args.group == "cases":
            covariates = covariates[group_values == str(args.cases_value)].copy()
        elif args.group == "controls":
            covariates = covariates[group_values != str(args.cases_value)].copy()

        if args.group_column in covariates.columns:
            covariates = covariates.drop(columns=[args.group_column])

    if args.group == "combined" and args.group_column:
        sample_group_column = resolve_column_name(samples.columns, args.group_column)
        covariate_group_column = resolve_column_name(covariates.columns, args.group_column)

        if sample_group_column and args.group_column not in covariates.columns:
            group_frame = samples[["IID", sample_group_column]].copy()
            if sample_group_column != args.group_column:
                group_frame = group_frame.rename(columns={sample_group_column: args.group_column})
            covariates = covariates.merge(group_frame, on="IID", how="left", sort=False)
        elif covariate_group_column:
            if covariate_group_column != args.group_column:
                covariates = covariates.rename(columns={covariate_group_column: args.group_column})
            pass
        else:
            raise ValueError(f"Group column not found in sample or covariate file: {args.group_column}")

    covariate_columns = [col for col in covariates.columns if col not in {"FID", "IID"}]
    for column in list(covariate_columns):
        series = covariates[column]
        non_missing = series.dropna()

        if non_missing.empty:
            print(f"WARNING: Dropping all-NA covariate: {column}")
            covariates = covariates.drop(columns=[column])
            continue

        unique_values = pd.unique(non_missing)
        unique_count = len(unique_values)
        if unique_count < 2:
            print(f"WARNING: Dropping invariant covariate: {column}")
            covariates = covariates.drop(columns=[column])
            continue

        if pd.api.types.is_numeric_dtype(series):
            if unique_count == 2:
                sorted_values = sorted(unique_values)
                print(f">>> Info: Converting binary covariate {column} ({sorted_values[0]}, {sorted_values[1]}) to 0/1")
                covariates[column] = binary_to_zero_one(series)
            else:
                covariates[column] = series.fillna(series.mean())
            continue

        if unique_count > (len(covariates) * 0.5):
            print(f"WARNING: Dropping covariate with too many levels: {column}")
            covariates = covariates.drop(columns=[column])
            continue

        mode_value = non_missing.mode().iloc[0]
        filled = series.fillna(mode_value)
        dummies = pd.get_dummies(filled, prefix=column, prefix_sep="", drop_first=True, dtype=int)
        covariates = pd.concat([covariates.drop(columns=[column]), dummies], axis=1)

    covariates = covariates[["FID", "IID"] + [col for col in covariates.columns if col not in {"FID", "IID"}]]
    covariates.to_csv(outdir / "covariates.cov", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
