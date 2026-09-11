#!/bin/bash
set -euo pipefail

input_files=""
protein_id=""
study=""
group=""
metal_binary="metal"
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --input_files) input_files="$2"; shift 2 ;;
    --protein_id) protein_id="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --metal_binary) metal_binary="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"
tmpdir=$(mktemp -d -p "$outdir")
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
metal_qc_script="${script_dir}/metal_qc.py"

if [[ ! -x "$metal_binary" ]]; then
  echo "ERROR: METAL binary not found or not executable: $metal_binary" > "${outdir}/meta.log"
  exit 1
fi

if [[ ! -f "$metal_qc_script" ]]; then
  echo "ERROR: METAL QC helper not found: $metal_qc_script" > "${outdir}/meta.log"
  exit 1
fi

IFS=',' read -r -a FILES <<< "$input_files"
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: No input files were provided to METAL" > "${outdir}/meta.log"
  exit 1
fi

input_list="${tmpdir}/metal_input.txt"
standardized_header=$'CHR\tPOS\tID\tEA\tOA\tEAF\tBETA\tSE\tP\tN'

read_header() {
  local infile="$1"

  if gzip -t "$infile" >/dev/null 2>&1; then
    python3 - "$infile" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", encoding="utf-8", errors="replace") as handle:
    print(handle.readline().rstrip("\r\n"))
PY
  else
    sed -n '1p' "$infile" | tr -d '\r'
  fi
}

# Stream the input (gz or plain) once and emit the METAL-ready file directly with
# the canonical position-based MARKER id (CHR:POS:sorted-alleles), while also
# appending the canonical-key -> original-rsID pair to rsidfile. This folds the
# former separate "canonicalize" rewrite and rsID-harvest passes into formatting.
# The canonical key makes cross-study matching unambiguous (rsIDs collide on
# multiallelics); the unsorted EA/OA columns still let METAL align effect
# direction, and CHR/POS are embedded so they are never lost.
reader_for() {
  if [[ "$2" == "1" ]]; then gzip -cd -- "$1"; else cat -- "$1"; fi
}

format_standardized_input() {
  local infile="$1" outfile="$2" compressed="$3" rsidfile="$4"
  reader_for "$infile" "$compressed" | awk -v rsidfile="$rsidfile" 'BEGIN{FS=OFS="\t"}
    NR==1 { print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"; next }
    NF == 0 { next }
    {
      chr=$1; gsub(/^[Cc][Hh][Rr]/, "", chr); pos=$2; oid=$3
      a=toupper($4); b=toupper($5)
      key=(a<=b)? chr":"pos":"a":"b : chr":"pos":"b":"a
      print key,$4,$5,$6,$7,$8,$9,$10,chr,pos
      if (oid ~ /^rs/ && chr!="" && pos!="") print key, oid >> rsidfile
    }' > "$outfile"
}

format_generic_input() {
  local infile="$1" outfile="$2" compressed="$3" rsidfile="$4"
  reader_for "$infile" "$compressed" | awk -v rsidfile="$rsidfile" 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF == 0) next
        if (!header_done) {
          for (i = 1; i <= NF; i++) {
            gsub(/"/, "", $i); gsub(/\r/, "", $i)
            t = toupper($i); gsub(/[^A-Z0-9]/, "", t)
            if (t == "ID" || t == "MARKERNAME" || t == "SNP") id = i
            else if (t == "EA" || t == "ALLELE1" || t == "A1" || t == "EFFECTALLELE") ea = i
            else if (t == "OA" || t == "ALLELE2" || t == "A2" || t == "OTHERALLELE") oa = i
            else if (t == "EAF" || t == "FREQ1" || t == "A1FREQ" || t == "FREQ" || t == "FREQA1") eaf = i
            else if (t == "BETA" || t == "EFFECT") beta = i
            else if (t == "SE" || t == "STDERR" || t == "STANDARDERROR") se = i
            else if (t == "P" || t == "PVALUE" || t == "PVAL") p = i
            else if (t == "N" || t == "WEIGHT" || t == "TOTALSAMPLESIZE" || t == "NEFF") n = i
            else if (t == "CHR" || t == "CHROM") chr = i
            else if (t == "POS" || t == "POSITION" || t == "GENPOS") pos = i
          }
          if (id == 0 || ea == 0 || oa == 0 || eaf == 0 || beta == 0 || se == 0 || p == 0) {
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
          header_done = 1
          next
        }
        nval = (n > 0 ? $(n) : "1")
        chrval = (chr > 0 ? $(chr) : "")
        posval = (pos > 0 ? $(pos) : "")
        gsub(/^[Cc][Hh][Rr]/, "", chrval)
        oid = $(id); aa = toupper($(ea)); bb = toupper($(oa))
        if (chrval != "" && posval != "")
          key = (aa <= bb) ? chrval":"posval":"aa":"bb : chrval":"posval":"bb":"aa
        else
          key = oid
        print key,$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval,chrval,posval
        if (oid ~ /^rs/ && chrval != "" && posval != "") print key, oid >> rsidfile
      }' > "$outfile"
}

input_data_rows=0
n_summary="${tmpdir}/n_summary.tsv"
printf 'ID\tN\n' > "$n_summary"
rsid_pairs="${tmpdir}/rsid_pairs.tsv"
: > "$rsid_pairs"
out_files=()

for i in "${!FILES[@]}"; do
  f="${FILES[$i]}"
  bn=$(basename "$f")
  out_txt="${tmpdir}/input_${i}_${bn%.gz}.txt"
  compressed="0"
  header=$(read_header "$f")

  if gzip -t "$f" >/dev/null 2>&1; then
    compressed="1"
  fi

  # Formatting now emits the canonical MARKER id and the rsID pairs in one pass.
  if [[ "$header" == "${standardized_header}"* ]]; then
    format_standardized_input "$f" "$out_txt" "$compressed" "$rsid_pairs"
  else
    format_generic_input "$f" "$out_txt" "$compressed" "$rsid_pairs"
  fi

  data_rows=$(( $(wc -l < "$out_txt") - 1 ))
  if (( data_rows <= 0 )); then
    echo "WARNING: METAL input has no variant rows after formatting: $f" >> "${outdir}/meta.log"
  fi
  input_data_rows=$((input_data_rows + data_rows))

  echo "$out_txt" >> "$input_list"
  out_files+=("$out_txt")
done

if (( input_data_rows == 0 )); then
  echo "ERROR: No variant rows were available for METAL after formatting inputs for ${protein_id} (${group})." >> "${outdir}/meta.log"
  exit 1
fi

# Per-marker summaries in a SINGLE pass over the formatted inputs:
#   n_summary     : summed sample size per id (METAL TotalSampleSize fallback)
#   coord_summary : authoritative CHR/POS per id (first occurrence; from REGENIE
#                   CHROM/GENPOS) -- METAL keeps only the MarkerName, so metal_qc
#                   joins CHR/POS back by id.
coord_summary="${tmpdir}/coord_summary.tsv"
awk -v nfile="$n_summary" -v cfile="$coord_summary" 'BEGIN{ FS=OFS="\t"
    print "ID","N" > nfile; print "ID","CHR","POS" > cfile }
  FNR==1 { next }
  NF > 0 && $1 != "" {
    n[$1] += $8
    if ($9 != "" && $10 != "" && !($1 in coord)) coord[$1] = $9 OFS $10
  }
  END {
    for (id in n)     print id, n[id]     > nfile
    for (id in coord) print id, coord[id] > cfile
  }' "${out_files[@]}"

# Canonical-key -> rsID map (first rsID seen per key), for re-attaching rsIDs to
# the meta output. Variants present only as chr:pos:ref:alt simply have no entry.
rsid_summary="${tmpdir}/rsid_summary.tsv"
printf 'ID\tRSID\n' > "$rsid_summary"
awk 'BEGIN { FS=OFS="\t" } $1 != "" && $2 != "" && !seen[$1]++ { print }' "$rsid_pairs" >> "$rsid_summary"

metal_script="${tmpdir}/metal_script.txt"
cat > "$metal_script" <<EOF
# TRACKPOSITIONS ON ensures METAL correctly tracks and outputs Chromosome/Position
TRACKPOSITIONS ON
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
MARKER ID
ALLELE EA OA
EFFECT BETA
STDERR SE
PVALUE P
WEIGHT N
FREQ EAF
CHROMOSOMELABEL CHR
POSITIONLABEL POS
EOF

while read -r f; do
  echo "PROCESS $f" >> "$metal_script"
done < "$input_list"

cat >> "$metal_script" <<EOF
OUTFILE ${outdir}/meta .tbl
ANALYZE HETEROGENEITY
QUIT
EOF

"$metal_binary" "$metal_script" > "${outdir}/meta.log" 2>&1

tbl_file=""
if [[ -f "${outdir}/meta.tbl" ]]; then
  tbl_file="${outdir}/meta.tbl"
else
  tbl_file=$(ls -1 "${outdir}"/meta*.tbl 2>/dev/null | head -n 1 || true)
fi

if [[ -z "$tbl_file" ]]; then
  echo "ERROR: METAL did not produce a .tbl output for ${protein_id}" >> "${outdir}/meta.log"
  exit 1
fi

stem="${tbl_file%.tbl}"
if [[ -f "${stem}.tbl.info" && "${stem}.tbl.info" != "${outdir}/meta.tbl.info" ]]; then
  mv -f "${stem}.tbl.info" "${outdir}/meta.tbl.info"
fi

python3 "$metal_qc_script" \
  --input_file "$tbl_file" \
  --protein_id "$protein_id" \
  --study "$study" \
  --group "$group" \
  --outdir "$outdir" \
  --n_file "$n_summary" \
  --coord_file "$coord_summary" \
  --rsid_file "$rsid_summary" \
  --require_heterogeneity

# Diagnostic: keep the raw METAL output header + a few rows so we can confirm
# which columns METAL emits (notably whether TRACKPOSITIONS produced
# Chromosome/Position) and what the MarkerName now looks like after re-keying.
head -n 6 "$tbl_file" > "${outdir}/meta.raw.head.tsv" 2>/dev/null || true

rm -f "$tbl_file"

echo "METAL for $protein_id ($group) complete."
