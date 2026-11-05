#!/bin/bash
#SBATCH --job-name=009_prep
#SBATCH --output=src/logs/009_prepare.log
#SBATCH --error=src/logs/009_prepare.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=high_p
#SBATCH --open-mode=truncate
#
# Assemble everything the permutation null needs, and PROVE the plumbing before
# any permutation runs. METHODS.md 8.3.
#
# Per study it collects:
#   covariates.cov   from the 'all' group  -- every sample, Sex/Age_Recr/PHENO/PC1-10
#   full.pheno       from the 'all' group  -- every sample, all proteins
#   keep_samples.txt from ever/ and never/ -- the TRUE stratum labels
#   dosages at the index variants          -- extracted by POSITION, not by ID
#
# Extraction is by position because the per-study .pvar and the summary-stat SNP
# key (chr:pos) need not use the same variant IDs.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${SLURM_SUBMIT_DIR:-${GWAS_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}}"
_LIB="${PROJ_ROOT}/src/lib/env.sh"
[[ -f "$_LIB" ]] || { echo "ERROR: ${_LIB} not found; sync the tree." >&2; exit 1; }
# shellcheck disable=SC1090
source "$_LIB"; load_env_file "${PROJ_ROOT}/.env" || true
mkdir -p "${PROJ_ROOT}/src/logs"

RESULTS=""; CMP=""; OUT="${PROJ_ROOT}/analysis/permutation"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results)    RESULTS="$2"; shift 2 ;;
    --comparison) CMP="$2"; shift 2 ;;
    --out)        OUT="$2"; shift 2 ;;
    -h|--help)    sed -n '10,24p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$RESULTS" ]] || { echo "ERROR: --results <pipeline outdir> required" >&2; exit 1; }
[[ -n "$CMP" ]] || CMP="${PROJ_ROOT}/analysis/comparison/pqtl_comparison.tsv"
[[ -f "$CMP" ]] || { echo "ERROR: not found: $CMP" >&2; exit 1; }
GEN="${GWAS_GENETICS_DIR:?GWAS_GENETICS_DIR not set}"
mkdir -p "$OUT/geno"

echo "=========================================="
echo " 009 permutation null -- input assembly"
echo "=========================================="
echo "Results   : $RESULTS"
echo "Comparison: $CMP"
echo "Genetics  : $GEN"
echo "Out       : $OUT"
echo

# ---- 1. the variant x protein pairs to refit ------------------------------
# PRIMARY, QC-passing loci only: those are the ones the global tests use.
python3 - "$CMP" "$OUT" <<'PY'
import csv, sys, pathlib
cmp_path, out = sys.argv[1], pathlib.Path(sys.argv[2])
pairs, variants = [], {}
with open(cmp_path) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r.get('TIER') != 'primary': continue
        if str(r.get('QC_PASS')).lower() not in ('true','1'): continue
        snp = r['SNP']
        try:
            chrom, pos = snp.split(':')[0], int(snp.split(':')[1])
        except (ValueError, IndexError):
            continue
        pairs.append((r['PROTEIN'], snp, r.get('EA',''), r.get('OA','')))
        variants[snp] = (chrom, pos)
with (out/'pairs.tsv').open('w') as fh:
    fh.write('protein\tsnp\tea\toa\n')
    for pr, s, ea, oa in pairs: fh.write(f'{pr}\t{s}\t{ea}\t{oa}\n')
# BED for plink2 --extract bed1 (0-based start, end exclusive)
with (out/'index_variants.bed').open('w') as fh:
    for snp,(c,p) in sorted(variants.items(), key=lambda kv:(kv[1][0], kv[1][1])):
        fh.write(f'{c}\t{p-1}\t{p}\t{snp}\n')
print(f"  pairs (protein x variant) : {len(pairs)}")
print(f"  distinct index variants   : {len(variants)}")
print(f"  distinct proteins         : {len(set(r[0] for r in pairs))}")
PY
echo

# ---- 2. per-study inputs ---------------------------------------------------
PLINK2="$(command -v plink2 || ls "${GWAS_ENV_DIR:-/nonexistent}"/*/bin/plink2 2>/dev/null | head -1)"
[[ -n "$PLINK2" ]] || { echo "ERROR: plink2 not found (tried PATH and \$GWAS_ENV_DIR/*/bin)" >&2; exit 1; }
echo "plink2: $PLINK2"; echo

printf '%-12s %8s %8s %8s %9s %9s %s\n' STUDY ALL EVER NEVER PROTEINS VARIANTS STATUS
ok=0; bad=0
for d in "$RESULTS"/*/; do
  study="$(basename "$d")"
  # Only real study directories have the per-group _shared tree. 'meta' and
  # 'pipeline_info' are pipeline output, not studies.
  [[ "$study" == _* ]] && continue
  [[ "$study" == "meta" || "$study" == "pipeline_info" ]] && continue
  [[ -d "$d/_shared" ]] || continue
  cov="$d/_shared/all/003_prepare-covariates/covariates.cov"
  phe="$d/_shared/all/002_prepare-phenotypes/full.pheno"
  kev="$d/_shared/ever/002_prepare-phenotypes/keep_samples.txt"
  knv="$d/_shared/never/002_prepare-phenotypes/keep_samples.txt"
  miss=""
  for f in "$cov" "$phe" "$kev" "$knv"; do [[ -f "$f" ]] || miss="$miss $(basename "$(dirname "$f")")/$(basename "$f")"; done
  if [[ -n "$miss" ]]; then
    printf '%-12s %8s %8s %8s %9s %9s MISSING:%s\n' "$study" - - - - - "$miss"; bad=$((bad+1)); continue
  fi
  # locate the fileset: flat <dir>/<study>.pgen or nested <dir>/<study>/<study>.pgen
  pfx=""
  [[ -f "$GEN/$study.pgen" ]] && pfx="$GEN/$study"
  [[ -f "$GEN/$study/$study.pgen" ]] && pfx="$GEN/$study/$study"
  if [[ -z "$pfx" ]]; then
    printf '%-12s %8s %8s %8s %9s %9s NO PGEN\n' "$study" - - - - -; bad=$((bad+1)); continue
  fi
  nall=$(( $(wc -l < "$cov") - 1 )); nev=$(wc -l < "$kev"); nnv=$(wc -l < "$knv")
  nprot=$(( $(head -1 "$phe" | tr '\t' '\n' | wc -l) - 2 ))
  # No --export-allele: it takes an ID/allele table, not a BED, so it was being
  # silently ignored. Export counts the ALT allele and the .raw header records
  # which allele that was; the companion .pvar gives REF/ALT per variant, and the
  # engine harmonises against the summary-stat EA/OA from pairs.tsv.
  "$PLINK2" --pfile "$pfx" \
      --extract bed1 "$OUT/index_variants.bed" \
      --export A \
      --out "$OUT/geno/$study" >/dev/null 2>"$OUT/geno/$study.plinkerr" || true
  "$PLINK2" --pfile "$pfx" \
      --extract bed1 "$OUT/index_variants.bed" \
      --make-just-pvar \
      --out "$OUT/geno/$study" >/dev/null 2>>"$OUT/geno/$study.plinkerr" || true
  # --- stage covariates, phenotypes and the TRUE stratum labels -----------
  mkdir -p "$OUT/cov" "$OUT/pheno" "$OUT/labels"
  cp -f "$cov" "$OUT/cov/$study.cov"
  cp -f "$phe" "$OUT/pheno/$study.pheno"
  # labels: IID, true stratum, and the permutation CELL. Within a study the cell
  # is case status, so shuffling inside it preserves both the study-level
  # stratum imbalance and the smoking-disease association (METHODS.md 8.3).
  python3 - "$cov" "$kev" "$knv" "$OUT/labels/$study.tsv" <<'PYL'
import sys
cov, kev, knv, out = sys.argv[1:5]
def ids(p):
    s = set()
    for i, ln in enumerate(open(p)):
        f = ln.split()
        if not f: continue
        if i == 0 and f[-1].upper() in ("IID", "ID"): continue
        s.add(f[-1])
    return s
ever, never = ids(kev), ids(knv)
head = open(cov).readline().split()
try:
    ip, ii = head.index("PHENO"), head.index("IID")
except ValueError:
    ip, ii = -1, 1
n = 0
with open(out, "w") as fh:
    fh.write("iid\tsmoke\tcell\n")
    for i, ln in enumerate(open(cov)):
        if i == 0: continue
        f = ln.split()
        if not f: continue
        iid = f[ii]
        arm = "EVER" if iid in ever else ("NEVER" if iid in never else None)
        if arm is None: continue          # unknown smoking: in neither arm
        cell = f[ip] if ip >= 0 else "0"
        fh.write(f"{iid}\t{arm}\t{cell}\n"); n += 1
print(f"      labels: {n} samples ({len(ever)} ever / {len(never)} never in keep lists)")
PYL

  if [[ -f "$OUT/geno/$study.raw" ]]; then
    nvar=$(( $(head -1 "$OUT/geno/$study.raw" | tr '\t' '\n' | wc -l) - 6 ))
    printf '%-12s %8s %8s %8s %9s %9s ok\n' "$study" "$nall" "$nev" "$nnv" "$nprot" "$nvar"; ok=$((ok+1))
  else
    printf '%-12s %8s %8s %8s %9s %9s PLINK FAILED (see geno/%s.plinkerr)\n' \
      "$study" "$nall" "$nev" "$nnv" "$nprot" - "$study"; bad=$((bad+1))
  fi
done
echo
for f in "$OUT"/geno/*.pvar; do
  [[ -f "$f" ]] || continue
  s="$(basename "${f%.pvar}")"
  dup=$(awk '!/^#/{print $1":"$2}' "$f" | sort | uniq -d | wc -l | tr -d ' ')
  [[ "$dup" -gt 0 ]] && echo "  note: $s has $dup multi-allelic position(s) -- engine harmonises on EA/OA"
done
echo
echo "studies ready: $ok   failed: $bad"
[[ $bad -eq 0 ]] || { echo "Fix the failures above before running the permutations." >&2; exit 1; }
echo "Inputs assembled under $OUT"
