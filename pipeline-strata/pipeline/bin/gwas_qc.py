#!/usr/bin/env python3
import argparse
import os
import numpy as np
import polars as pl
from scipy.stats import chi2

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_files", required=True, help="Comma-separated REGENIE files")
    parser.add_argument("--protein_id", required=True)
    parser.add_argument("--info_score", type=float, default=0.3)
    parser.add_argument("--study", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()

def calculate_lambda(p_values):
    # Filter out nulls
    p = p_values.drop_nulls()
    if len(p) == 0: return np.nan
    # Clip to avoid math errors (1e-300)
    p_arr = np.clip(p.to_numpy(), 1e-300, 1.0)
    return np.nanmedian(chi2.isf(p_arr, df=1)) / chi2.ppf(0.5, df=1)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    files = [f.strip() for f in args.input_files.split(",")]
    
    # Use Polars to read and combine
    dfs = []
    for f in files:
        try:
            df = pl.read_csv(f, separator=" ")
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")
            continue
    
    if not dfs:
        print("Error: No input files could be read.")
        exit(1)
        
    df = pl.concat(dfs)
    
    # Immediate subset to save memory
    keep_cols = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "LOG10P", "N", "BETA", "SE"]
    df = df.select([c for c in keep_cols if c in df.columns])

    # Harmonize and Calculate P
    df = df.with_columns([
        (10.0 ** -pl.col("LOG10P")).alias("P")
    ])
    
    df = df.rename({
        "CHROM": "CHR",
        "GENPOS": "POS",
        "ALLELE0": "OA",
        "ALLELE1": "EA",
        "A1FREQ": "EAF"
    })

    l_raw = calculate_lambda(df["P"])

    # QC Filter (INFO already filtered upstream at info_min; keep N >= 50% max)
    max_n = df["N"].max()
    if max_n is not None and max_n > 0:
        df_filt = df.filter(pl.col("N") >= 0.5 * max_n)
    else:
        df_filt = df
        
    l_filt = calculate_lambda(df_filt["P"])

    # Strictly enforce the output column order as requested by user
    # Order: [CHR, POS, ID, EA, OA, EAF, BETA, SE, P, N, INFO]
    standard_cols = ["CHR", "POS", "ID", "EA", "OA", "EAF", "BETA", "SE", "P", "N", "INFO"]
    
    # Ensure CHR is stripped of "chr" prefix safely (cast to string first)
    if "CHR" in df_filt.columns:
        df_filt = df_filt.with_columns(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", ""))


    final_df = df_filt.select([c for c in standard_cols if c in df_filt.columns])
    
    # High-speed sink_csv() streaming.
    gwas_out = os.path.join(args.outdir, "gwas.tsv.gz")
    final_df.lazy().sink_csv(gwas_out, separator="\t", compression="gzip")


    # Metrics
    metrics_path = os.path.join(args.outdir, "metrics.tsv")
    metrics_df = pl.DataFrame([{
        "study": args.study,
        "group": args.group,
        "id": args.protein_id, 
        "l_raw": l_raw, 
        "l_filt": l_filt
    }])
    metrics_df.write_csv(metrics_path, separator="\t")

if __name__ == "__main__":
    main()
