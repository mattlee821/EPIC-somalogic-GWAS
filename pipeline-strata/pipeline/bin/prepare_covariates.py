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
    # See prepare_phenotypes.py: stratum membership is an explicit whitelist,
    # never the complement of the case value.
    parser.add_argument("--group_values", default="")
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

    # Capture the stratifying column BEFORE covariate subsetting, which would
    # otherwise drop it. It may live in either the sample file or the covariate
    # file; the sample file wins so a .psam column is never shadowed.
    group_frame = None
    if args.group_column.strip():
        sample_group_column = resolve_column_name(samples.columns, args.group_column)
        covariate_group_column = resolve_column_name(covariates.columns, args.group_column)

        if sample_group_column:
            group_frame = samples[["IID", sample_group_column]].copy()
            group_frame.columns = ["IID", args.group_column]
        elif covariate_group_column:
            group_frame = covariates[["IID", covariate_group_column]].copy()
            group_frame.columns = ["IID", args.group_column]
        else:
            raise ValueError(
                f"Group column '{args.group_column}' not found in the sample file or the "
                f"covariate file. Sample columns: {sorted(samples.columns)[:20]}. "
                f"Covariate columns: {sorted(covariates.columns)[:20]}."
            )

    if args.include_covariates.strip():
        requested_covariates = [part.strip() for part in args.include_covariates.split(",") if part.strip()]

        # A requested covariate may live in the .psam rather than the covariate
        # file. The frozen pipeline only ever pulled the GROUP column from the
        # .psam, so disease status (PHENO, which exists ONLY in the .psam --
        # 003 deliberately keeps it out of the phenofile) was reachable as a
        # covariate only while it was also the grouping variable. Stratified
        # runs group on smoking/BMI and must still adjust for PHENO, so any
        # requested covariate now falls back to the sample file.
        for name in requested_covariates:
            if resolve_column_name(covariates.columns, name):
                continue
            psam_column = resolve_column_name(samples.columns, name)
            if psam_column is None:
                continue
            extra = samples[["IID", psam_column]].copy()
            extra.columns = ["IID", name]
            covariates = covariates.merge(extra, on="IID", how="left", sort=False)
            print(f">>> Info: covariate {name} resolved from the sample file (.psam)")

        # Fail loudly rather than silently fitting a different model than asked
        # for. A study whose .psam carries no PHENO column will trip this; drop
        # PHENO from --covariates or drop the study, but make it a decision.
        missing = [c for c in requested_covariates if resolve_column_name(covariates.columns, c) is None]
        if missing:
            raise ValueError(
                f"Requested covariates not found in the covariate file or the sample file: "
                f"{missing}. Available covariate columns: {sorted(covariates.columns)[:30]}. "
                f"Available sample columns: {sorted(samples.columns)[:20]}."
            )

        resolved = [resolve_column_name(covariates.columns, c) for c in requested_covariates]
        covariates = covariates[["FID", "IID"] + resolved].copy()
        covariates.columns = ["FID", "IID"] + requested_covariates
    else:
        exclude_columns = {"STUDY", "FID", "IID", "FID_1", "IID_1", "PlateId"}
        proteins_to_exclude = read_path_list(args.include_proteins)
        if proteins_to_exclude is not None:
            exclude_columns.update(proteins_to_exclude)
        # Never let the stratifying variable leak in as a covariate: within a
        # stratum it is constant, and in the unstratified reference group it
        # would be a different model than the strata are fitted with.
        if args.group_column.strip():
            exclude_columns.add(args.group_column)

        numeric_columns = [col for col in covariates.columns if pd.api.types.is_numeric_dtype(covariates[col])]
        keep_columns = [col for col in numeric_columns if col not in exclude_columns]
        covariates = covariates[["FID", "IID"] + keep_columns].copy()

    if group_frame is not None:
        accepted = [v.strip() for v in args.group_values.split(",") if v.strip()]
        if not accepted:
            raise ValueError(
                f"Group '{args.group}' declares group_column '{args.group_column}' but no "
                "--group_values. See the 'strata' block in params.yaml."
            )

        covariates = covariates.merge(group_frame, on="IID", how="inner", sort=False)
        observed = covariates[args.group_column].astype("string")
        n_before = len(covariates)
        covariates = covariates[observed.isin(accepted)].copy()
        print(
            f">>> Stratum {args.group}: {args.group_column} in {accepted} -> "
            f"{len(covariates)}/{n_before} samples"
        )
        covariates = covariates.drop(columns=[args.group_column])

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

    # Record the model actually fitted. A covariate can be dropped here for one
    # study and kept for another -- PHENO is invariant in the all-case studies
    # (Ecvd_01, Ecvd_02) and survives elsewhere -- which means a meta-analysis
    # can silently pool studies fitted on different models. Writing the surviving
    # covariates next to the results makes that auditable afterwards instead of
    # being visible only in a Nextflow work directory.
    final_covariates = [col for col in covariates.columns if col not in {"FID", "IID"}]
    requested = [c.strip() for c in args.include_covariates.split(",") if c.strip()]
    dropped = [c for c in requested if resolve_column_name(final_covariates, c) is None]
    with (outdir / "covariates_used.txt").open("w") as handle:
        handle.write(f"study\t{args.study_id}\n")
        handle.write(f"group\t{args.group}\n")
        handle.write(f"n_samples\t{len(covariates)}\n")
        handle.write(f"requested\t{','.join(requested) if requested else '<auto>'}\n")
        handle.write(f"fitted\t{','.join(final_covariates)}\n")
        handle.write(f"dropped\t{','.join(dropped) if dropped else '<none>'}\n")
    if dropped:
        print(f">>> NOTE: requested covariate(s) dropped for {args.study_id}/{args.group}: {dropped}")


if __name__ == "__main__":
    main()
