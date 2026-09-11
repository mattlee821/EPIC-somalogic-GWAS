#!/bin/bash
set -euo pipefail

POSITIONAL_ARGS=()
pfile=""
study=""
king=""
keep=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pfile="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --king) king="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$pfile" || -z "$study" || -z "$king" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <prefix> --study <val> --king <val> --outdir <val> [--keep <file>]" >&2
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
log_file="${outdir}/sample_qc.log"
exec > >(tee -a "$log_file") 2>&1

echo "=========================================="
echo " Starting sample QC for $study"
echo " Input PLINK2 prefix: $pfile"
echo "=========================================="

KEEP_ARG=""
if [[ -n "$keep" && -f "$keep" ]]; then
    echo "Restricting analysis to samples in: $keep"
    KEEP_ARG="--keep $keep"
fi

# Step 1: Heterozygosity
echo ">>> Step 1: Heterozygosity check"
plink2 --pfile "$pfile" \
  $KEEP_ARG \
  --indep-pairwise 200 50 0.25 \
  --out "${outdir}/het_ld_temp"

plink2 --pfile "$pfile" \
  $KEEP_ARG \
  --extract "${outdir}/het_ld_temp.prune.in" \
  --het \
  --out "${outdir}/het"

# Use Python to avoid an extra R/data.table dependency in the PLINK env
echo "Filtering heterozygosity outliers (mean +/- 3SD)..."
python - <<EOF
import csv
import math
from pathlib import Path

het_path = Path("${outdir}/het.het")
out_path = Path("${outdir}/het_pass.txt")

with het_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)

if not rows:
    raise SystemExit("No rows found in het.het")

fieldnames = rows[0].keys()
f_col = "F" if "F" in fieldnames else list(fieldnames)[-1]
iid_col = next((col for col in fieldnames if col.upper().endswith("IID")), None)
fid_col = next((col for col in fieldnames if col.upper().endswith("FID")), None)

if iid_col is None:
    raise SystemExit("Could not find IID column in het.het")

f_vals = []
for row in rows:
    value = row.get(f_col, "")
    if value not in ("", "NA", "nan", "NaN"):
        f_vals.append(float(value))

if not f_vals:
    raise SystemExit("No valid F values found in het.het")

mean_f = sum(f_vals) / len(f_vals)
if len(f_vals) > 1:
    variance = sum((value - mean_f) ** 2 for value in f_vals) / (len(f_vals) - 1)
    sd_f = math.sqrt(variance)
else:
    sd_f = 0.0

lower = mean_f - 3 * sd_f
upper = mean_f + 3 * sd_f

print(f"F-coeff mean: {mean_f:.4f}, SD: {sd_f:.4f}")
print(f"Keep range: [{lower:.4f}, {upper:.4f}]")

with out_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    for row in rows:
        value = row.get(f_col, "")
        if value in ("", "NA", "nan", "NaN"):
            continue
        f_value = float(value)
        if lower <= f_value <= upper:
            if fid_col:
                writer.writerow([row[fid_col], row[iid_col]])
            else:
                writer.writerow([row[iid_col]])
EOF

# Step 2: Relatedness filter (--king-cutoff $king)
echo ">>> Step 2: Relatedness filter (--king-cutoff $king)"
plink2 --pfile "$pfile" \
  --keep "${outdir}/het_pass.txt" \
  --king-cutoff "$king" \
  --out "${outdir}/king_pass"

# Final clean sample list
if [[ -f "${outdir}/king_pass.king.cutoff.in.id" ]]; then
    cp "${outdir}/king_pass.king.cutoff.in.id" "${outdir}/qc_pass_samples.txt"
else
    # In case no relatives found to prune
    cp "${outdir}/het_pass.txt" "${outdir}/qc_pass_samples.txt"
fi

final_count=$(wc -l < "${outdir}/qc_pass_samples.txt")
echo "=========================================="
echo " QC Finished: $final_count samples passed all filters."
echo " Final list: ${outdir}/qc_pass_samples.txt"
echo "=========================================="
