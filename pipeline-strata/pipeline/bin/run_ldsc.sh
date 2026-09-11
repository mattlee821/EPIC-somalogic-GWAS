#!/bin/bash
set -euo pipefail

# run_ldsc.sh -- per-protein LD Score Regression intercept + SNP-h2 on one
# cross-study meta sumstats file (meta.tsv.gz from META_STUDY). Emits a single
# tab-separated row (with header) to --out_row so the Nextflow LDSC module can
# collectFile() them into metrics_ldsc.tsv.
#
# The intercept distinguishes genuine polygenic signal / power (intercept ~1)
# from confounding or stratification (intercept >> 1); SNP-h2 (Observed scale)
# is reported as a by-product -- this is the same "polygenic component" estimator
# the UKB-PPP paper used (LDSC, Bulik-Sullivan 2015), not GREML.
#
# The meta sumstats are keyed by chr:pos:sorted-alleles but carry an RSID column;
# LDSC needs rsIDs, so we use RSID directly and restrict to HapMap3. We do NOT
# pass --merge-alleles: the intercept depends only on chi2 and LD score
# (allele-orientation invariant), so the HM3 SNP list is a sufficient filter.
#
# One protein never fails the pipeline: on any munge/ldsc error we still write an
# NA row and exit 0, so the protein appears in metrics_ldsc.tsv.

meta=""
protein_id=""
ldsc_home=""
ref_dir=""
hm3=""
outdir="."
out_row=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --meta) meta="$2"; shift 2 ;;
    --protein_id) protein_id="$2"; shift 2 ;;
    --ldsc_home) ldsc_home="$2"; shift 2 ;;
    --ref_dir) ref_dir="$2"; shift 2 ;;
    --hm3) hm3="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --out_row) out_row="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$outdir"
[[ -n "$out_row" ]] || out_row="${outdir}/ldsc.tsv"
log="${outdir}/ldsc.run.log"

header=$'protein\tn_snps\tintercept\tintercept_se\th2\th2_se\tlambda_gc'

# Emit an NA row and exit cleanly; used for any non-fatal per-protein failure so
# the protein still appears in the consolidated metrics_ldsc.tsv.
emit_na() {
  local reason="$1" nsnp="${2:-NA}"
  echo "WARNING: ${protein_id}: ${reason}" | tee -a "$log" >&2
  printf '%s\n' "$header" > "$out_row"
  printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\n' "$protein_id" "$nsnp" >> "$out_row"
  exit 0
}

for p in "${ldsc_home}/munge_sumstats.py" "${ldsc_home}/ldsc.py" "$hm3"; do
  [[ -e "$p" ]] || emit_na "missing LDSC asset $p (run: bash src/004_ldsc-herit.sh --reference)"
done
[[ -d "$ref_dir" ]] || emit_na "LD ref dir missing: $ref_dir (run --setup)"
[[ -f "$meta" ]] || emit_na "meta sumstats not found: $meta"

# Confirm the meta carries an RSID column (read header once; avoid grep -q SIGPIPE
# tripping pipefail).
meta_header="$(zcat "$meta" 2>/dev/null | head -1)" || true
case $'\t'"${meta_header}"$'\t' in
  *$'\t'RSID$'\t'*) : ;;
  *) emit_na "meta has no RSID column; cannot map to HapMap3 rsIDs" ;;
esac

prep="${outdir}/munge_input.txt"
# Build munge input (SNP A1 A2 N Z P), restricted to HapMap3, rsID from RSID.
# Z = BETA/SE (signed). Skip rows with empty RSID (e.g. indels) or SE<=0.
zcat "$meta" | awk -v hm3f="$hm3" '
  BEGIN{ FS=OFS="\t"
    while ((getline s < hm3f) > 0) { if (s != "SNP") hm3[s]=1 } }
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; print "SNP","A1","A2","N","Z","P"; next }
  {
    snp = $(h["RSID"])
    if (snp == "" || !(snp in hm3)) next
    se = $(h["SE"]) + 0
    if (se > 0) print snp, $(h["EA"]), $(h["OA"]), $(h["N"]), $(h["BETA"])/se, $(h["P"])
  }' > "$prep"
nsnp=$(( $(wc -l < "$prep") - 1 ))
echo "${protein_id}: HapMap3 SNPs for munge: ${nsnp}" | tee -a "$log"
(( nsnp >= 200 )) || emit_na "too few HapMap3 SNPs (${nsnp}) for a reliable intercept" "$nsnp"

# Munge (no --merge-alleles: intercept is allele-orientation invariant).
python "${ldsc_home}/munge_sumstats.py" \
  --sumstats "$prep" \
  --snp SNP --a1 A1 --a2 A2 --N-col N --p P --signed-sumstats Z,0 \
  --chunksize 500000 \
  --out "${outdir}/munged" >> "$log" 2>&1 \
  || emit_na "munge_sumstats failed (see ${log})" "$nsnp"

# LDSC --h2 (intercept is the quantity of interest; h2 reported as by-product).
python "${ldsc_home}/ldsc.py" \
  --h2 "${outdir}/munged.sumstats.gz" \
  --ref-ld-chr "${ref_dir}/" \
  --w-ld-chr "${ref_dir}/" \
  --out "${outdir}/h2" >> "$log" 2>&1 \
  || emit_na "ldsc --h2 failed (see ${outdir}/h2.log)" "$nsnp"

hlog="${outdir}/h2.log"
icpt_line=$(grep -iE "^Intercept:" "$hlog" 2>/dev/null | tail -1 || true)
h2_line=$(grep -iE "Total Observed scale h2:|Observed scale h2:" "$hlog" 2>/dev/null | tail -1 || true)
lam_line=$(grep -iE "Lambda GC:" "$hlog" 2>/dev/null | tail -1 || true)
icpt=$(echo "$icpt_line"   | sed -E 's/.*Intercept:[[:space:]]*([0-9.eE+-]+).*/\1/')
icpt_se=$(echo "$icpt_line" | sed -E 's/.*\(([0-9.eE+-]+)\).*/\1/')
h2=$(echo "$h2_line"       | sed -E 's/.*h2:[[:space:]]*([0-9.eE+-]+).*/\1/')
h2_se=$(echo "$h2_line"    | sed -E 's/.*\(([0-9.eE+-]+)\).*/\1/')
lam=$(echo "$lam_line"     | sed -E 's/.*Lambda GC:[[:space:]]*([0-9.eE+-]+).*/\1/')

printf '%s\n' "$header" > "$out_row"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$protein_id" "${nsnp}" "${icpt:-NA}" "${icpt_se:-NA}" "${h2:-NA}" "${h2_se:-NA}" "${lam:-NA}" >> "$out_row"
echo "${protein_id}: intercept=${icpt:-NA} (${icpt_se:-NA})  h2=${h2:-NA} (${h2_se:-NA})  lambda=${lam:-NA}" | tee -a "$log"
