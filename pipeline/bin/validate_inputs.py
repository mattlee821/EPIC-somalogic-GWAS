#!/usr/bin/env python3
import argparse
import sys
import os

from config import read_sample_file, resolve_column_name

def check_file_exists(filepath, desc):
    if not os.path.exists(filepath):
        print(f"FAIL: {desc} '{filepath}' does not exist (Current Dir: {os.getcwd()})")
        return False
    return True

def normalize_chromosome(value):
    cleaned = str(value).strip()
    if cleaned.lower().startswith("chr"):
        cleaned = cleaned[3:]
    return cleaned.upper()

def read_pvar_chromosomes(pvar_path):
    chroms = set()
    chrom_idx = 0
    header_seen = False

    with open(pvar_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue
            fields = line.split()
            if not header_seen and fields:
                maybe_header = [field.lstrip("#").upper() for field in fields]
                if "CHROM" in maybe_header:
                    chrom_idx = maybe_header.index("CHROM")
                    header_seen = True
                    continue
                header_seen = True
            if len(fields) <= chrom_idx:
                continue
            chroms.add(normalize_chromosome(fields[chrom_idx]))

    return chroms

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--phenotype_file", required=True)
    parser.add_argument("--covariate_file", required=True)
    parser.add_argument("--include_studies", default="")
    parser.add_argument("--include_proteins", default="")
    parser.add_argument("--covariates", default="")
    parser.add_argument("--chromosomes", default="")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    report_path = os.path.join(args.outdir, "validation_report.txt")
    
    error_found = False
    with open(report_path, "w") as f:
        def report(message, *, is_error=False):
            print(message, file=f)
            print(message, file=sys.stderr if is_error else sys.stdout)

        # 1. Samplesheet
        if not check_file_exists(args.samplesheet, "Samplesheet"):
            report("FAIL: Samplesheet missing", is_error=True)
            sys.exit(1)
            
        studies = []
        with open(args.samplesheet, "r") as s:
            header = s.readline().strip().split(",")
            required_columns = ["study_id", "pfile", "sample_file", "group_column", "cases_value"]
            missing_columns = [col for col in required_columns if col not in header]
            if missing_columns:
                report(
                    f"FAIL: Samplesheet missing required column(s): {', '.join(missing_columns)}",
                    is_error=True,
                )
                error_found = True
            if not missing_columns:
                for line in s:
                    row = line.strip().split(",")
                    if len(row) > 0:
                        studies.append(dict(zip(header, row)))
                    
        req_studies = [x.strip() for x in args.include_studies.split(",") if x.strip()] if args.include_studies else []
        requested_chromosomes = [
            normalize_chromosome(x)
            for x in args.chromosomes.split(",")
            if x.strip()
        ]

        samplesheet_studies = [s.get("study_id", "") for s in studies]
        if req_studies:
            missing_requested = [sid for sid in req_studies if sid not in samplesheet_studies]
            if missing_requested:
                report(
                    f"FAIL: Requested study/studies absent from samplesheet: {', '.join(missing_requested)}",
                    is_error=True,
                )
                error_found = True
            selected_studies = [s for s in studies if s.get("study_id") in req_studies]
        else:
            selected_studies = studies

        if not selected_studies:
            report("FAIL: No studies selected from samplesheet.", is_error=True)
            error_found = True
        
        # 2. Check study files
        for s in selected_studies:
            sid = s['study_id']
                
            report(f"Checking study {sid}...")
            
            pfile_prefix = s.get('pfile', '').strip()
            if not pfile_prefix:
                report(f"FAIL: Missing pfile value for {sid}", is_error=True)
                error_found = True
                continue

            missing_pfile_parts = [
                f"{pfile_prefix}.{ext}"
                for ext in ["pgen", "pvar", "psam"]
                if not os.path.exists(f"{pfile_prefix}.{ext}")
            ]
            if missing_pfile_parts:
                report(
                    f"FAIL: Incomplete PLINK2 input for {sid}: missing {', '.join(missing_pfile_parts)}",
                    is_error=True,
                )
                error_found = True
            else:
                report(f"Found complete PLINK2 pfile for {sid}")

            sample_file = s.get('sample_file', '').strip()
            expected_sample_file = os.path.abspath(f"{pfile_prefix}.psam")
            observed_sample_file = os.path.abspath(sample_file) if sample_file else ""
            if not sample_file:
                report(f"FAIL: Missing sample_file value for {sid}", is_error=True)
                error_found = True
            elif observed_sample_file != expected_sample_file:
                report(
                    f"FAIL: sample_file for {sid} must equal <pfile>.psam. Expected {expected_sample_file}, observed {observed_sample_file}",
                    is_error=True,
                )
                error_found = True
            elif not os.path.exists(sample_file):
                report(f"FAIL: sample_file for {sid} does not exist: {sample_file}", is_error=True)
                error_found = True

            if requested_chromosomes and os.path.exists(f"{pfile_prefix}.pvar"):
                observed_chromosomes = read_pvar_chromosomes(f"{pfile_prefix}.pvar")
                missing_chromosomes = [
                    chrom for chrom in requested_chromosomes
                    if chrom not in observed_chromosomes
                ]
                # A missing autosome is a real data problem (hard fail). A missing
                # sex/mito chromosome (X/23, Y/24, MT/26) is tolerated: some cohorts
                # legitimately lack it (e.g. sex-specific cancer studies) and
                # REGENIE step 2 guards per study, skipping absent chromosomes.
                sex_mito = {"X", "23", "Y", "24", "MT", "M", "26"}
                missing_autosomes = [c for c in missing_chromosomes if c not in sex_mito]
                missing_sexmito = [c for c in missing_chromosomes if c in sex_mito]
                if missing_autosomes:
                    report(
                        f"FAIL: Requested autosome(s) absent from {pfile_prefix}.pvar for {sid}: {', '.join(missing_autosomes)}",
                        is_error=True,
                    )
                    error_found = True
                if missing_sexmito:
                    report(
                        f"WARNING: Requested chromosome(s) absent from {pfile_prefix}.pvar for {sid}: {', '.join(missing_sexmito)} (skipped for this study)",
                        is_error=False,
                    )

        # 3. Phenotype Header
        if not check_file_exists(args.phenotype_file, "Phenotype file"):
            report("FAIL: Missing phenotype file", is_error=True)
            error_found = True
        else:
            with open(args.phenotype_file, "r") as pf:
                line = pf.readline().strip()
                pheno_header = line.split("\t") if "\t" in line else line.split(",")
                
                req_proteins = []
                if args.include_proteins:
                    if os.path.exists(args.include_proteins):
                        with open(args.include_proteins, "r") as f_prots:
                            req_proteins = [x.strip() for x in f_prots if x.strip()]
                    else:
                        req_proteins = [x.strip() for x in args.include_proteins.split(",")]
                
                if req_proteins:
                    missing_proteins = [p for p in req_proteins if p not in pheno_header]
                    if missing_proteins:
                        # Only report first 10 missing to avoid flooding
                        show_missing = missing_proteins[:10]
                        report(
                            f"FAIL: Missing requested proteins (showing {len(show_missing)}/{len(missing_proteins)}): {', '.join(show_missing)}",
                            is_error=True,
                        )
                        error_found = True

        # 4. Covariate Header (Fallback to phenotype file if covariates.txt doesn't exist)
        cov_source = args.covariate_file if os.path.exists(args.covariate_file) else args.phenotype_file
        
        if os.path.exists(cov_source):
            with open(cov_source, "r") as cf:
                line = cf.readline().strip()
                cov_header = line.split("\t") if "\t" in line else line.split(",")
                
                if args.covariates:
                    req_covs = [x.strip() for x in args.covariates.split(",")]
                    missing_covs = [c for c in req_covs if resolve_column_name(cov_header, c) is None]
                    if missing_covs:
                        report(f"FAIL: Missing requested covariates: {', '.join(missing_covs)}", is_error=True)
                        error_found = True
        else:
            report(f"FAIL: Covariate source '{cov_source}' not found.", is_error=True)
            error_found = True
            cov_header = []

        for s in selected_studies:
            # Only check group_column if it's actually specified in samplesheet
            gcol = s.get('group_column', '').strip()
            if not gcol:
                continue

            in_cov = resolve_column_name(cov_header, gcol) is not None
            in_sample = False
            sample_file = s.get('sample_file', '')
            if sample_file and os.path.exists(sample_file):
                try:
                    sample_df = read_sample_file(sample_file)
                except Exception as exc:
                    report(f"FAIL: Could not read sample file for {s['study_id']}: {exc}", is_error=True)
                    error_found = True
                else:
                    if resolve_column_name(sample_df.columns, gcol) is not None:
                        in_sample = True

            if not in_cov and not in_sample:
                report(
                    f"FAIL: Group column '{gcol}' missing for {s['study_id']} (covariate/sample files)",
                    is_error=True,
                )
                error_found = True
        
        if error_found:
            report("Validation finished with ERRORS.", is_error=True)
        else:
            report("Validation PASS.")

    if error_found:
        print("Validation finished with ERRORS. Check reports.")
        sys.exit(1)
    else:
        print("Validation PASS.")

if __name__ == "__main__":
    main()
