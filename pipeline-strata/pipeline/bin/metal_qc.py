#!/usr/bin/env python3
import argparse
import os
import tempfile
import numpy as np
import polars as pl
from scipy.stats import chi2

COLUMN_MAPPING = {
    "MarkerName": "ID",
    "Allele1": "EA",
    "Allele2": "OA",
    "Freq1": "EAF",
    "Effect": "BETA",
    "StdErr": "SE",
    "P-value": "P",
    "TotalSampleSize": "N",
    "Weight": "N",
    "N": "N",
    "Chromosome": "CHR",
    "Position": "POS",
    "CHROMOSOME": "CHR",
    "POSITION": "POS",
    "Direction": "DIRECTION",
    "HetISq": "HET_ISQ",
    "HetChiSq": "HET_CHISQ",
    "HetDf": "HET_DF",
    "HetPVal": "HET_P",
}

STANDARD_COLS = ["CHR", "POS", "ID", "RSID", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
META_COLS = ["DIRECTION", "HET_ISQ", "HET_CHISQ", "HET_DF", "HET_P"]

def normalize_col_name(name):
    return "".join(ch for ch in name.upper() if ch.isalnum())

NORMALIZED_COLUMN_MAPPING = {
    normalize_col_name(old): new
    for old, new in COLUMN_MAPPING.items()
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True, help="Single METAL output file")
    parser.add_argument("--protein_id", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--n_file", default="")
    parser.add_argument("--coord_file", default="",
                        help="TSV (ID, CHR, POS) of authoritative coordinates "
                             "from the per-study inputs, joined by ID.")
    parser.add_argument("--rsid_file", default="",
                        help="TSV (ID, RSID) mapping the canonical chr:pos:alleles "
                             "key back to the original rsID, joined by ID.")
    parser.add_argument("--require_heterogeneity", action="store_true")
    return parser.parse_args()

def calculate_lambda(p_values):
    p = p_values.drop_nulls()
    if len(p) == 0: return np.nan
    p_arr = np.clip(p.to_numpy(), 1e-300, 1.0)
    return np.nanmedian(chi2.isf(p_arr, df=1)) / chi2.ppf(0.5, df=1)

def has_expected_columns(columns):
    known = set(COLUMN_MAPPING) | set(COLUMN_MAPPING.values()) | set(STANDARD_COLS) | set(META_COLS)
    normalized_known = set(NORMALIZED_COLUMN_MAPPING) | {normalize_col_name(c) for c in known}
    return bool(set(columns) & known) or bool({normalize_col_name(c) for c in columns} & normalized_known)

def read_metal_table(input_file):
    """Read METAL output whether it is tab-delimited or whitespace-delimited."""
    df = pl.read_csv(input_file, separator="\t", ignore_errors=True, infer_schema_length=10000)
    if has_expected_columns(df.columns):
        return df

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            delete=False,
            dir=os.path.dirname(os.path.abspath(input_file)),
        ) as tmp:
            tmp_path = tmp.name
            with open(input_file, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    stripped = line.strip()
                    if stripped:
                        tmp.write("\t".join(stripped.split()))
                        tmp.write("\n")

        return pl.read_csv(tmp_path, separator="\t", ignore_errors=True, infer_schema_length=10000)
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)

def rename_known_columns(df):
    rename_map = {}
    existing_targets = set(df.columns)
    scheduled_targets = set()

    for old in df.columns:
        new = COLUMN_MAPPING.get(old)
        if new is None:
            new = NORMALIZED_COLUMN_MAPPING.get(normalize_col_name(old))
        if new is None or old == new:
            continue
        if new in existing_targets or new in scheduled_targets:
            continue
        rename_map[old] = new
        scheduled_targets.add(new)

    if rename_map:
        df = df.rename(rename_map)
    return df

def require_columns(df, columns, context):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: {context} missing required columns after parsing: {','.join(missing)}"
        )

    all_null = [c for c in columns if df[c].null_count() == df.height]
    if all_null:
        raise SystemExit(
            f"ERROR: {context} columns are present but entirely empty: {','.join(all_null)}"
        )

def fill_n_from_file(df, n_file):
    if not n_file:
        return df
    if not os.path.exists(n_file):
        raise SystemExit(f"ERROR: N summary file not found: {n_file}")

    n_df = pl.read_csv(n_file, separator="\t", ignore_errors=True)
    require_columns(n_df, ["ID", "N"], "N summary file")
    if n_df.height == 0:
        raise SystemExit("ERROR: N summary file contains no variant rows.")
    n_df = n_df.select([
        pl.col("ID").cast(pl.Utf8),
        pl.col("N").cast(pl.Float64, strict=False).alias("N_FROM_INPUTS"),
    ])

    df = df.with_columns(pl.col("ID").cast(pl.Utf8))
    df = df.join(n_df, on="ID", how="left")
    if "N" in df.columns:
        df = df.with_columns(
            pl.coalesce([
                pl.col("N").cast(pl.Float64, strict=False),
                pl.col("N_FROM_INPUTS"),
            ]).alias("N")
        )
    else:
        df = df.rename({"N_FROM_INPUTS": "N"})

    if "N_FROM_INPUTS" in df.columns:
        df = df.drop("N_FROM_INPUTS")
    return df

def fill_coords_from_file(df, coord_file):
    """Join authoritative CHR/POS by ID from the per-study inputs.

    METAL keeps only the marker ID, from which positions cannot be recovered for
    rsIDs (and only fragilely for chr:pos:ref:alt IDs). The per-study sumstats
    carry reliable CHR/POS (REGENIE CHROM/GENPOS), so we join those back by ID and
    let them take precedence over anything parsed from the ID string.
    """
    if not coord_file or not os.path.exists(coord_file):
        return df

    c_df = pl.read_csv(coord_file, separator="\t", ignore_errors=True)
    if not {"ID", "CHR", "POS"}.issubset(c_df.columns):
        return df

    c_df = c_df.select([
        pl.col("ID").cast(pl.Utf8),
        pl.col("CHR").cast(pl.Utf8).str.replace("(?i)^chr", "").alias("CHR_REF"),
        pl.col("POS").cast(pl.Int64, strict=False).alias("POS_REF"),
    ]).unique(subset=["ID"], keep="first")

    df = df.with_columns(pl.col("ID").cast(pl.Utf8)).join(c_df, on="ID", how="left")

    chr_sources = ["CHR_REF"] + (["CHR"] if "CHR" in df.columns else [])
    pos_sources = ["POS_REF"] + (["POS"] if "POS" in df.columns else [])
    df = df.with_columns([
        pl.coalesce([pl.col(c).cast(pl.Utf8) for c in chr_sources]).alias("CHR"),
        pl.coalesce([pl.col(c).cast(pl.Int64, strict=False) for c in pos_sources]).alias("POS"),
    ]).drop(["CHR_REF", "POS_REF"])
    return df

def fill_rsid_from_file(df, rsid_file):
    """Re-attach the original rsID by joining on the canonical chr:pos:alleles ID.

    The meta MarkerName is the positional key (rsIDs collide on multiallelics, so
    they cannot be the matching key); this restores the rsID as a separate column
    for annotation/lookup. Variants with no rsID (chr:pos:ref:alt only) stay null.
    """
    if not rsid_file or not os.path.exists(rsid_file):
        return df
    r_df = pl.read_csv(rsid_file, separator="\t", ignore_errors=True)
    if not {"ID", "RSID"}.issubset(r_df.columns):
        return df
    r_df = r_df.select([
        pl.col("ID").cast(pl.Utf8),
        pl.col("RSID").cast(pl.Utf8),
    ]).unique(subset=["ID"], keep="first")
    df = df.with_columns(pl.col("ID").cast(pl.Utf8)).join(r_df, on="ID", how="left")
    return df

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Read the raw METAL output or an already-standardized table, then keep only
    # the columns needed for final output and QC metrics.
    df = read_metal_table(args.input_file)
    df = rename_known_columns(df)
    keep_cols = STANDARD_COLS + META_COLS
    df = df.select([c for c in keep_cols if c in df.columns])

    required_cols = ["ID", "EA", "OA", "EAF", "BETA", "SE", "P"]
    require_columns(df, required_cols, "METAL output")

    if df.height == 0:
        raise SystemExit("ERROR: METAL output contains no variant rows.")

    if args.require_heterogeneity:
        require_columns(df, META_COLS, "METAL heterogeneity output")

    df = fill_n_from_file(df, args.n_file)
    df = fill_coords_from_file(df, args.coord_file)
    df = fill_rsid_from_file(df, args.rsid_file)

    cast_exprs = []
    for col in ("P", "N", "EAF", "BETA", "SE", "HET_ISQ", "HET_CHISQ", "HET_DF", "HET_P"):
        if col in df.columns:
            cast_exprs.append(pl.col(col).cast(pl.Float64, strict=False))
    if "POS" in df.columns:
        cast_exprs.append(pl.col("POS").cast(pl.Int64, strict=False))
    if cast_exprs:
        df = df.with_columns(cast_exprs)

    if args.require_heterogeneity:
        require_columns(df, META_COLS, "METAL heterogeneity output")

    # Ensure CHR and POS columns are present and correctly handled
    needs_split = False
    if "CHR" not in df.columns or "POS" not in df.columns:
        needs_split = True
    elif df["CHR"].null_count() == len(df) or df["POS"].null_count() == len(df):
        needs_split = True

    if needs_split and "ID" in df.columns:
        # replace_all (not replace, which only swaps the first ':') so a
        # chr:pos:ref:alt ID splits into [chr, pos, ref, alt] and POS = field 1.
        id_clean = df["ID"].cast(pl.Utf8).str.replace_all(":", "_")
        id_parts = id_clean.str.split("_")

        new_cols = []
        if "CHR" not in df.columns:
            new_cols.append(id_parts.list.get(0, null_on_oob=True).str.replace("(?i)chr", "").alias("CHR"))
        elif df["CHR"].null_count() > 0:
            # Fill nulls only
            new_cols.append(pl.when(pl.col("CHR").is_null())
                            .then(id_parts.list.get(0, null_on_oob=True).str.replace("(?i)chr", ""))
                            .otherwise(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", ""))
                            .alias("CHR"))
        else:
             new_cols.append(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", "").alias("CHR"))

        if "POS" not in df.columns:
            new_cols.append(id_parts.list.get(1, null_on_oob=True).cast(pl.Int64, strict=False).alias("POS"))
        elif df["POS"].null_count() > 0:
            new_cols.append(pl.when(pl.col("POS").is_null())
                            .then(id_parts.list.get(1, null_on_oob=True).cast(pl.Int64, strict=False))
                            .otherwise(pl.col("POS"))
                            .alias("POS"))
        
        df = df.with_columns(new_cols)
    else:
        if "CHR" in df.columns:
            df = df.with_columns(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", ""))

    if "EA" in df.columns:
        df = df.with_columns(pl.col("EA").cast(pl.Utf8).str.to_uppercase())
    if "OA" in df.columns:
        df = df.with_columns(pl.col("OA").cast(pl.Utf8).str.to_uppercase())

    l_raw = calculate_lambda(df["P"])

    if "N" in df.columns:
        max_n = df["N"].max()
        if max_n is not None and max_n > 0:
            df_filt = df.filter(pl.col("N") >= 0.5 * max_n)
        else:
            df_filt = df
    else:
        df_filt = df
    
    l_filt = calculate_lambda(df_filt["P"])

    output_cols = STANDARD_COLS + META_COLS

    for col in output_cols:
        if col not in df_filt.columns:
            df_filt = df_filt.with_columns(pl.lit(None).alias(col))

    final_df = df_filt.select(output_cols)

    if final_df.height == 0:
        raise SystemExit("ERROR: No variant rows remain after meta-analysis QC filtering.")

    gwas_out = os.path.join(args.outdir, "meta.tsv.gz")
    final_df.lazy().sink_csv(gwas_out, separator="\t", compression="gzip")

    metrics_path = os.path.join(args.outdir, "metrics.tsv")
    metrics_df = pl.DataFrame([{
        "study": args.study,
        "group": args.group,
        "id": args.protein_id,
        "n_raw": df.height,
        "n_filt": final_df.height,
        "l_raw": l_raw,
        "l_filt": l_filt
    }])
    metrics_df.write_csv(metrics_path, separator="\t")

if __name__ == "__main__":
    main()
