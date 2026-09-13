#!/bin/bash
set -euo pipefail

# run_herit.sh -- per-protein heritability decomposition on one cross-study meta
# sumstats file, following UKB-PPP (Sun et al. 2023):
#   total SNP-h2 = pQTL component + polygenic component
#     pQTL component   : sum of variance explained (2*f*(1-f)*beta^2) by the
#                        independent index pQTLs from clumping, split cis/trans.
#     polygenic compt. : LDSC --h2 over the genome EXCLUDING the pQTL regions
#                        (index SNP +/- window) and the HLA region.
# Emits a single tab-separated row (with header) to --out_row for collectFile().
#
# Robust: any failure writes an NA-filled row and exits 0 so the protein still
# appears in metrics_herit.tsv and never kills the pipeline.

meta="" protein_id="" ld_ref="" protein_map=""
ldsc_home="" ref_dir="" hm3=""
cis_window_kb="1000" clump_p="5e-8" clump_r2="0.1"
hla="6:25500000:34000000"
outdir="." out_row=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --meta) meta="$2"; shift 2 ;;
    --protein_id) protein_id="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --ld_ref) ld_ref="$2"; shift 2 ;;
    --protein_map) protein_map="$2"; shift 2 ;;
    --ldsc_home) ldsc_home="$2"; shift 2 ;;
    --ref_dir) ref_dir="$2"; shift 2 ;;
    --hm3) hm3="$2"; shift 2 ;;
    --cis_window_kb) cis_window_kb="$2"; shift 2 ;;
    --clump_p) clump_p="$2"; shift 2 ;;
    --clump_r2) clump_r2="$2"; shift 2 ;;
    --hla) hla="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --out_row) out_row="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$outdir"
[[ -n "$out_row" ]] || out_row="${outdir}/herit.tsv"
log="${outdir}/herit.run.log"
# See run_ldsc.sh: the stratum belongs in the row, or the join cross-products.
header=$'protein\tgroup\tn_index\tn_cis\tn_trans\th2_cis\th2_trans\th2_polygenic\th2_total\tintercept\tintercept_se\tlambda_gc'

emit_na() {
  echo "WARNING: ${protein_id}: $1" | tee -a "$log" >&2
  printf '%s\n' "$header" > "$out_row"
  printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$protein_id" "${group:-NA}" >> "$out_row"
  exit 0
}

[[ -f "$meta" ]] || emit_na "meta sumstats not found: $meta"
[[ -f "$protein_map" ]] || emit_na "protein_gene_map not found: $protein_map (run: bash src/004_ldsc-herit.sh --protein-map)"

# Header check: meta must carry the columns we use.
meta_header="$(zcat "$meta" 2>/dev/null | head -1)" || true
for need in CHR POS RSID EA OA EAF BETA SE P N; do
  case $'\t'"${meta_header}"$'\t' in *$'\t'"$need"$'\t'*) : ;; *) emit_na "meta missing column $need" ;; esac
done

# MHC/HLA handling: the MHC's long-range, atypical LD breaks both clumping
# (spurious "independent" index pQTLs) and the LD-Score model, so -- matching the
# LDSC convention (we munge against hm3_no_MHC) and UKB-PPP (HLA excluded from
# clumping) -- we drop the HLA window from EVERY component here. Reported h2 is
# therefore MHC-excluded and consistent across the pQTL and polygenic parts.
IFS=':' read -r hla_c hla_a hla_b <<< "$hla"

# 1. Deduped assoc table (one row per chr:pos, lowest P) for clumping + variance,
#    with the MHC/HLA region removed so it cannot inflate the pQTL component.
assoc="${outdir}/assoc.tsv"
zcat "$meta" | awk -v hc="$hla_c" -v ha="$hla_a" -v hb="$hla_b" 'BEGIN{FS=OFS="\t"}
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  {
    chr=$(h["CHR"]); gsub(/^[Cc][Hh][Rr]/,"",chr); pos=$(h["POS"])
    if (chr=="" || pos=="") next
    if (chr==hc && pos+0>=ha+0 && pos+0<=hb+0) next   # exclude MHC/HLA
    snp=chr":"pos; p=$(h["P"])+0
    if (!(snp in seenP) || p < seenP[snp]) {
      seenP[snp]=p
      row[snp]=snp OFS chr OFS pos OFS $(h["EA"]) OFS $(h["OA"]) OFS $(h["EAF"]) OFS $(h["BETA"]) OFS p
    }
  }
  END{ print "SNP","CHR","POS","EA","OA","EAF","BETA","P"; for(s in row) print row[s] }' > "$assoc"

# 2. Independent index pQTLs. Prefer LD clumping against the combined reference;
#    fall back to distance-based greedy pruning if no reference is available.
index="${outdir}/index_snps.txt"
window_kb_int=${cis_window_kb%.*}
if [[ -n "$ld_ref" && -f "${ld_ref}.bed" ]]; then
  echo "${protein_id}: clumping against LD reference ${ld_ref}" | tee -a "$log"
  if plink --bfile "$ld_ref" \
       --clump "$assoc" --clump-snp-field SNP --clump-field P \
       --clump-p1 "$clump_p" --clump-p2 "$clump_p" \
       --clump-r2 "$clump_r2" --clump-kb "$window_kb_int" \
       --out "${outdir}/clumped" >> "$log" 2>&1 && [[ -f "${outdir}/clumped.clumped" ]]; then
    awk 'NR>1 && $3!="" {print $3}' "${outdir}/clumped.clumped" > "$index"
  else
    echo "${protein_id}: plink --clump produced no output (no genome-wide hits?)" | tee -a "$log"
    : > "$index"
  fi
else
  echo "${protein_id}: no LD reference; distance-based greedy pruning (window ${window_kb_int}kb)" | tee -a "$log"
  # Greedy: sort by P, take a SNP as index if no kept index within +/-window on same chr.
  awk -v w=$((window_kb_int*1000)) 'BEGIN{FS=OFS="\t"}
    NR==1{next}
    { print $8, $2, $3, $1 }' "$assoc" \
    | sort -k1,1g \
    | awk -v w=$((window_kb_int*1000)) -v p="$clump_p" '{
        if (($1+0) > (p+0)) next
        chr=$2; pos=$3; snp=$4; keep=1
        for (i in kc) if (kc[i]==chr && (pos-kp[i] <= w && kp[i]-pos <= w)) { keep=0; break }
        if (keep){ n++; kc[n]=chr; kp[n]=pos; print snp }
      }' > "$index"
fi
n_index=$(wc -l < "$index" | tr -d ' ')
echo "${protein_id}: ${n_index} independent index pQTLs" | tee -a "$log"

# 3. Variance-explained (cis/trans) + the regions to exclude from the polygenic fit.
summary="${outdir}/pqtl_summary.txt"
exclude="${outdir}/exclude_regions.bed"
python "${script_dir}/herit_compute.py" \
  --assoc "$assoc" --index "$index" --protein-map "$protein_map" \
  --protein-id "$protein_id" --cis-window-kb "$cis_window_kb" \
  --out-summary "$summary" --out-exclude "$exclude" >> "$log" 2>&1 \
  || emit_na "herit_compute failed (see ${log})"
# shellcheck disable=SC1090
source "$summary"   # n_index, n_cis, n_trans, h2_cis, h2_trans, gene_chr
# Belt-and-suspenders: also drop the HLA window from the polygenic remainder
# (hm3_no_MHC already removes it; hla_c/a/b were parsed above).
printf '%s\t%s\t%s\n' "$hla_c" "$hla_a" "$hla_b" >> "$exclude"

# 4. Polygenic component: LDSC --h2 on HapMap3 SNPs OUTSIDE every excluded region.
poly_h2="NA"; intercept="NA"; intercept_se="NA"; lambda="NA"
if [[ -n "$ldsc_home" && -f "${ldsc_home}/munge_sumstats.py" && -d "$ref_dir" && -f "$hm3" ]]; then
  prep="${outdir}/munge_input.txt"
  zcat "$meta" | awk -v hm3f="$hm3" -v exf="$exclude" '
    BEGIN{ FS=OFS="\t"
      while ((getline s < hm3f) > 0) if (s!="SNP") hm3[s]=1
      ne=0; while ((getline e < exf) > 0){ split(e,a,"\t"); ne++; ec[ne]=a[1]; ea[ne]=a[2]+0; eb[ne]=a[3]+0 } }
    NR==1{ for(i=1;i<=NF;i++) h[$i]=i; print "SNP","A1","A2","N","Z","P"; next }
    {
      snp=$(h["RSID"]); if (snp=="" || !(snp in hm3)) next
      chr=$(h["CHR"]); gsub(/^[Cc][Hh][Rr]/,"",chr); pos=$(h["POS"])+0
      for(i=1;i<=ne;i++) if (chr==ec[i] && pos>=ea[i] && pos<=eb[i]) next
      se=$(h["SE"])+0; if (se>0) print snp,$(h["EA"]),$(h["OA"]),$(h["N"]),$(h["BETA"])/se,$(h["P"])
    }' > "$prep"
  nsnp=$(( $(wc -l < "$prep") - 1 ))
  echo "${protein_id}: polygenic-remainder HapMap3 SNPs: ${nsnp}" | tee -a "$log"
  if (( nsnp >= 200 )) \
     && python "${ldsc_home}/munge_sumstats.py" --sumstats "$prep" \
          --snp SNP --a1 A1 --a2 A2 --N-col N --p P --signed-sumstats Z,0 \
          --chunksize 500000 --out "${outdir}/munged" >> "$log" 2>&1 \
     && python "${ldsc_home}/ldsc.py" --h2 "${outdir}/munged.sumstats.gz" \
          --ref-ld-chr "${ref_dir}/" --w-ld-chr "${ref_dir}/" \
          --out "${outdir}/h2" >> "$log" 2>&1; then
    hlog="${outdir}/h2.log"
    il=$(grep -iE "^Intercept:" "$hlog" | tail -1 || true)
    hl=$(grep -iE "Total Observed scale h2:|Observed scale h2:" "$hlog" | tail -1 || true)
    ll=$(grep -iE "Lambda GC:" "$hlog" | tail -1 || true)
    poly_h2=$(echo "$hl"     | sed -E 's/.*h2:[[:space:]]*([0-9.eE+-]+).*/\1/'); poly_h2="${poly_h2:-NA}"
    intercept=$(echo "$il"   | sed -E 's/.*Intercept:[[:space:]]*([0-9.eE+-]+).*/\1/'); intercept="${intercept:-NA}"
    intercept_se=$(echo "$il"| sed -E 's/.*\(([0-9.eE+-]+)\).*/\1/'); intercept_se="${intercept_se:-NA}"
    lambda=$(echo "$ll"      | sed -E 's/.*Lambda GC:[[:space:]]*([0-9.eE+-]+).*/\1/'); lambda="${lambda:-NA}"
  else
    echo "${protein_id}: polygenic LDSC step failed/insufficient SNPs" | tee -a "$log"
  fi
else
  echo "${protein_id}: LDSC assets missing; polygenic component = NA" | tee -a "$log"
fi

# 5. Total = cis + trans + polygenic (treat a missing polygenic as 0 only if the
#    pQTL component succeeded; if everything is NA the row stays NA).
total=$(awk -v c="$h2_cis" -v t="$h2_trans" -v p="$poly_h2" 'BEGIN{
  s=c+t; if (p!="NA") s+=p; printf "%.6g", s }')

printf '%s\n' "$header" > "$out_row"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$protein_id" "${group:-NA}" "${n_index:-0}" "${n_cis:-0}" "${n_trans:-0}" \
  "${h2_cis:-NA}" "${h2_trans:-NA}" "${poly_h2}" "${total}" \
  "${intercept}" "${intercept_se}" "${lambda}" >> "$out_row"
echo "${protein_id}: cis=${h2_cis} trans=${h2_trans} poly=${poly_h2} total=${total}" | tee -a "$log"
