#!/bin/bash
set -euo pipefail

pfile=""
study=""
maf=""
info_min=""
outdir=""

# Sample exclusions (ancestry/related) and HWE flagging are performed upstream by
# the imputation pipeline and already applied to the merged pfile by 002. Variant
# QC here therefore only applies the MAF and imputation-quality (R2) filters.
while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pfile="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --maf) maf="$2"; shift 2 ;;
    --info_min) info_min="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$pfile" || -z "$study" || -z "$maf" || -z "$info_min" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <prefix> --study <val> --maf <val> --info_min <val> --outdir <val>"
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
log_file="${outdir}/variant_qc.log"
exec > >(tee -a "$log_file") 2>&1

echo "Starting variant QC for $study using PLINK2 prefix $pfile..."

plink2 --pfile "$pfile" \
  --maf "$maf" \
  --extract-if-info "R2 >= $info_min" \
  --write-snplist \
  --out "${outdir}/variant_qc_pass"

echo "Variant QC complete. Variant list written to ${outdir}/variant_qc_pass.snplist"
