#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import re

import pandas as pd

from config import ensure_string_ids, read_path_list, read_sample_file, read_table_auto, resolve_column_name, standardize_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--study_id", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--group_column", default="")
    parser.add_argument("--cases_value", default="")
    parser.add_argument("--covariate_file", required=True)
    parser.add_argument("--covariates", default="")
    parser.add_argument("--include_proteins", default="")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def build_feature_stem(index: int, protein_id: str) -> str:
    safe_token = re.sub(r"[^A-Za-z0-9._-]+", "_", protein_id.strip())
    safe_token = safe_token.strip("._-") or "trait"
    return f"feature_{index:05d}_{safe_token}"


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pheno = standardize_ids(read_table_auto(args.phenotype_file))
    samples = read_sample_file(args.sample_file)

    if "IID" not in pheno.columns or "IID" not in samples.columns:
        raise ValueError("Both phenotype and sample inputs must contain an IID column.")

    pheno = ensure_string_ids(pheno, ("IID",))
    samples = ensure_string_ids(samples, ("IID",))

    # Restrict to the samples assigned to THIS study by the cross-study dedup
    # (003). The phenotype file lists each person once with their assigned STUDY;
    # a person genotyped in several studies shares the same IID (Idepic_Bio) across
    # those studies' .psam files, so without this filter they would match every
    # such study by IID and be analysed — and counted in the meta N — multiple
    # times, silently undoing the dedup.
    if "STUDY" in pheno.columns:
        pheno = pheno[pheno["STUDY"].astype(str) == str(args.study_id)].copy()
        if pheno.empty:
            raise ValueError(f"No phenotype rows assigned to study '{args.study_id}'.")

    common_samples = pd.Index(pheno["IID"].dropna().unique()).intersection(samples["IID"].dropna().unique())
    if common_samples.empty:
        raise ValueError("Zero overlap found.")

    pheno = pheno[pheno["IID"].isin(common_samples)].copy()
    # The PLINK2 .psam from the imputation pipeline uses FID == IID, and plink2
    # --keep / REGENIE match on FID+IID. Setting FID = IID (not "0") makes the
    # keep, phenotype, and covariate files join correctly against the genotypes.
    pheno.insert(0, "FID", pheno["IID"])
    pheno = pheno[["FID", "IID"] + [col for col in pheno.columns if col not in {"FID", "IID"}]]

    all_columns = list(pheno.columns)
    id_columns = ["FID", "IID"]
    exclude_columns = {"STUDY", "FID", "IID", "FID_1", "IID_1"}

    if args.covariates.strip():
        exclude_columns.update(part.strip() for part in args.covariates.split(",") if part.strip())

    requested_proteins = read_path_list(args.include_proteins)

    if requested_proteins is not None:
        valid_proteins = [protein for protein in requested_proteins if protein in all_columns and protein not in exclude_columns]
        pheno = pheno[id_columns + valid_proteins].copy()
    else:
        numeric_columns = [col for col in pheno.columns if pd.api.types.is_numeric_dtype(pheno[col])]
        trait_columns = [col for col in numeric_columns if col not in exclude_columns]
        pheno = pheno[id_columns + trait_columns].copy()

    if args.group and args.group != "combined":
        covariates = standardize_ids(read_table_auto(args.covariate_file))
        if "IID" not in covariates.columns:
            raise ValueError("Covariate file must contain an IID column.")
        covariates = ensure_string_ids(covariates, ("IID",))

        sample_group_column = resolve_column_name(samples.columns, args.group_column)
        covariate_group_column = resolve_column_name(covariates.columns, args.group_column)

        if sample_group_column:
            group_frame = samples[["IID", sample_group_column]].copy()
            if sample_group_column != args.group_column:
                group_frame = group_frame.rename(columns={sample_group_column: args.group_column})
        elif covariate_group_column:
            group_frame = covariates[["IID", covariate_group_column]].copy()
            if covariate_group_column != args.group_column:
                group_frame = group_frame.rename(columns={covariate_group_column: args.group_column})
        else:
            raise ValueError(f"Group column not found in sample or covariate file: {args.group_column}")

        pheno = pheno.merge(group_frame, on="IID", how="inner", sort=False)
        group_values = pheno[args.group_column].astype("string")

        if args.group == "cases":
            pheno = pheno[group_values == str(args.cases_value)].copy()
        elif args.group == "controls":
            pheno = pheno[group_values != str(args.cases_value)].copy()

        pheno = pheno.drop(columns=[args.group_column])
        pheno = pheno[["FID", "IID"] + [col for col in pheno.columns if col not in {"FID", "IID"}]]

    full_out = outdir / "full.pheno"
    pheno.to_csv(full_out, sep="\t", index=False, na_rep="NA")

    keep_out = outdir / "keep_samples.txt"
    pheno[["FID", "IID"]].to_csv(keep_out, sep="\t", index=False, header=False)

    protein_columns = [col for col in pheno.columns if col not in {"FID", "IID"}]
    if not protein_columns:
        raise ValueError("No phenotype traits remain after filtering.")

    summary_rows = []
    manifest_rows = []

    for index, protein in enumerate(protein_columns, start=1):
        feature_out = outdir / f"{build_feature_stem(index, protein)}.pheno"
        pheno[id_columns + [protein]].to_csv(feature_out, sep="\t", index=False, na_rep="NA")
        manifest_rows.append({"protein_id": protein, "pheno_file": feature_out.name})

        values = pd.to_numeric(pheno[protein], errors="coerce")
        summary_rows.append(
            {
                "protein": protein,
                "n_total": len(values),
                "n_missing": int(values.isna().sum()),
                "mean": values.mean(),
                "sd": values.std(ddof=1),
                "min": values.min(),
                "max": values.max(),
            }
        )

    manifest_df = pd.DataFrame(manifest_rows, columns=["protein_id", "pheno_file"])
    manifest_df.to_csv(outdir / "features.manifest", sep="\t", index=False)

    summary_df = pd.DataFrame(
        summary_rows,
        columns=["protein", "n_total", "n_missing", "mean", "sd", "min", "max"],
    )
    summary_df.to_csv(outdir / "protein_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
