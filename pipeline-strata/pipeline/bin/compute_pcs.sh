#!/bin/bash
set -euo pipefail

step1_pfile=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) step1_pfile="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$step1_pfile" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <step1_input> --outdir <outdir>"
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${step1_pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${step1_pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"

plink2 --pfile "$step1_pfile" \
  --pca 10 \
  --out "${outdir}/pca"

if [[ ! -f "${outdir}/pca.eigenvec" ]]; then
  echo "ERROR: PCA failed to produce pca.eigenvec" >&2
  exit 1
fi
