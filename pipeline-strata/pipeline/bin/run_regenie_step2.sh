#!/bin/bash
set -euo pipefail

pgen_prefix=""
sample_file=""
pheno_file=""
cov_file=""
pred_list=""
study=""
group=""
chromosomes=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pgen_prefix="$2"; shift 2 ;;
    --sample_file) sample_file="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --pred_list) pred_list="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --chr|--chromosomes) chromosomes="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"

if [[ -z "$pgen_prefix" || -z "$sample_file" || -z "$pheno_file" || -z "$cov_file" || -z "$pred_list" || -z "$study" || -z "$group" || -z "$bsize" || -z "$outdir" ]]; then
    echo "Usage: $0 --pfile <prefix> --sample_file <psam> --pheno_file <file> --cov_file <file> --pred_list <file> --study <id> --group <group> --chromosomes <list> --bsize <val> --outdir <dir>" >&2
    exit 1
fi

if [[ -z "$chromosomes" ]]; then
    echo "ERROR: At least one chromosome must be provided via --chromosomes"
    exit 1
fi

for ext in pgen pvar psam; do
    if [[ ! -f "${pgen_prefix}.${ext}" ]]; then
        echo "ERROR: Missing PLINK2 input file: ${pgen_prefix}.${ext}" >&2
        exit 1
    fi
done

expected_sample_file="${pgen_prefix}.psam"
if [[ "$sample_file" != "$expected_sample_file" ]]; then
    echo "ERROR: sample_file must match the pfile .psam path." >&2
    echo "Expected: $expected_sample_file" >&2
    echo "Observed: $sample_file" >&2
    exit 1
fi

# Default to Nextflow-provided CPU count when available (e.g., NXF_TASK_CPUS).
# Fall back to SLURM_CPUS_PER_TASK, then 1 for manual runs.
threads="${NXF_TASK_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="${threads}"
export OPENBLAS_NUM_THREADS="${threads}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_MAX_THREADS="${threads}"

if [[ -n "$pred_list" ]]; then
    localized_pred_list="${outdir}/pred.list"
    python3 - "$pred_list" "$localized_pred_list" <<'PY'
import os
import re
import sys

src_path, dst_path = sys.argv[1:]
tmp_path = dst_path + ".tmp"
pattern = re.compile(r'\S+\.loco\.gz')

def rewrite_token(match: re.Match[str]) -> str:
    token = match.group(0)
    base = os.path.basename(token)
    if base.endswith(".loco.gz") and os.path.exists(base):
        return f"./{base}"
    return token

with open(src_path, "r", encoding="utf-8") as src, open(tmp_path, "w", encoding="utf-8") as dst:
    for raw_line in src:
        dst.write(pattern.sub(rewrite_token, raw_line))

os.replace(tmp_path, dst_path)
PY
    pred_list="${localized_pred_list}"
fi

# REGENIE 4.x requires #FID as the first column and IID as the second column.
ORIG_PSAM="${pgen_prefix}.psam"
# Read the first non-## line in a single awk (exit after match) rather than
# `grep -v '^##' | head -1`: with pipefail, head closing the pipe early sends
# SIGPIPE to grep on large psams, which would abort the script (exit 141).
header=$(awk '!/^##/{print; exit}' "$ORIG_PSAM")

if [[ "$header" != \#FID$'\t'IID* && "$header" != "#FID IID"* ]]; then
    echo ">>> Fixing PSAM header for REGENIE compatibility..."

    FIXED_PREFIX="${outdir}/fixed_gen"

    python3 - "$ORIG_PSAM" "${FIXED_PREFIX}.psam" <<'PY'
import sys

src_path, dst_path = sys.argv[1:]

with open(src_path, "r", encoding="utf-8") as src:
    raw_lines = [line.rstrip("\n") for line in src if line.strip()]

comment_lines = [line for line in raw_lines if line.startswith("##")]
data_lines = [line for line in raw_lines if not line.startswith("##")]

if not data_lines:
    raise SystemExit(f"No PSAM header found in {src_path}")

header = data_lines[0].split()
header[0] = header[0].lstrip("#")
upper = [col.upper() for col in header]

if "IID" not in upper:
    raise SystemExit(f"PSAM file {src_path} does not contain an IID column")

iid_idx = upper.index("IID")
fid_idx = upper.index("FID") if "FID" in upper else None
remaining_indices = [
    idx for idx in range(len(header))
    if idx != iid_idx and idx != fid_idx
]

out_header = ["#FID", "IID"] + [header[idx] for idx in remaining_indices]

with open(dst_path, "w", encoding="utf-8") as dst:
    for line in comment_lines:
        dst.write(line + "\n")
    dst.write("\t".join(out_header) + "\n")
    for line in data_lines[1:]:
        fields = line.split()
        if len(fields) < len(header):
            raise SystemExit(f"Malformed PSAM row in {src_path}: {line}")
        fid = fields[fid_idx] if fid_idx is not None else "0"
        iid = fields[iid_idx]
        out_fields = [fid, iid] + [fields[idx] for idx in remaining_indices]
        dst.write("\t".join(out_fields) + "\n")
PY

    ln -sf "$(readlink -f "${pgen_prefix}.pgen")" "${FIXED_PREFIX}.pgen"
    ln -sf "$(readlink -f "${pgen_prefix}.pvar")" "${FIXED_PREFIX}.pvar"

    pgen_prefix="${FIXED_PREFIX}"
fi

# Run REGENIE step 2 ONCE over the whole pgen, not once per chromosome. regenie
# processes every chromosome present in the pgen in a single pass, so a study that
# lacks a chromosome (e.g. X in some cohorts) simply has no variants for it and is
# handled without error -- no per-chromosome loop or guard is needed. The previous
# loop re-read the entire ~11M-variant pgen once per chromosome (~22x redundant
# NFS reads per task); under many concurrent step-2 jobs that saturated scratch IO
# and blew past the wall-time limit. One pass writes a single all-chromosome
# output (full_<phenotype>.regenie.gz), which GWAS_QC consumes directly.
# (--chromosomes is still accepted for interface compatibility but not used.)
log_file="${outdir}/step2.log"
cmd=(regenie --step 2 --threads "${threads}"
    --pgen "${pgen_prefix}"
    --phenoFile "${pheno_file}"
    --covarFile "${cov_file}"
    --pred "${pred_list}"
    --bsize "${bsize}"
    --qt
    --gz
    --out "${outdir}/full")

printf 'Running: %q ' "${cmd[@]}" | tee "$log_file"
printf '\n' | tee -a "$log_file"

if ! "${cmd[@]}" >> "$log_file" 2>&1; then
    echo "ERROR: REGENIE Step 2 failed. Log follows:"
    cat "$log_file"
    exit 1
fi

echo "REGENIE Step 2 for $study ($group) complete (all chromosomes in one pass)."
