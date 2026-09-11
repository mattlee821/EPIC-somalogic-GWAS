process GWAS_QC {
  label        'low'
  conda        "${params.env_dir}/py_qc"
  publishDir   params.outdir, mode: 'rellink', saveAs: { name ->
    def parts = name.tokenize('/')
    def protein_id = parts ? parts[0] : "UNKNOWN"
    def base = parts ? parts[-1] : name
    "${study_id}/${protein_id}/${group}/009_QC/${base}"
  }

  input:
  tuple val(study_id), val(group), val(protein_ids), path(regenie_files), val(mapping)

  output:
  tuple val(study_id), val(group), val(protein_ids), path("*/gwas.tsv.gz"), path("*/metrics.tsv"), emit: all_out
  path "*/gwas.tsv.gz", emit: gwas_files

  path "*/metrics.tsv", emit: metrics_files

  script:
  """
  export POLARS_MAX_THREADS="1"
  export RAYON_NUM_THREADS="1"
  export OMP_NUM_THREADS="1"
  export OPENBLAS_NUM_THREADS="1"
  export MKL_NUM_THREADS="1"
  export NUMEXPR_MAX_THREADS="1"
  parallel_jobs=${task.cpus}
  echo -e "${mapping}" > mapping.txt

  run_gwas_qc() {
    local prot="\$1"
    local files="\$2"
    mkdir -p "\${prot}"
    python ${projectDir}/bin/gwas_qc.py \\
      --input_files "\${files}" \\
      --protein_id "\${prot}" \\
      --study "${study_id}" \\
      --group "${group}" \\
      --info_score ${params.info_score} \\
      --outdir "\${prot}"
  }

  count_pids() {
    if [[ -z "\${pids}" ]]; then
      echo 0
    else
      printf '%s\\n' "\${pids}" | wc -l
    fi
  }

  reap_finished_pids() {
    local next_pids=""
    local pid

    while IFS= read -r pid; do
      [[ -n "\${pid}" ]] || continue
      if kill -0 "\${pid}" 2>/dev/null; then
        if [[ -n "\${next_pids}" ]]; then
          printf -v next_pids '%s\n%s' "\${next_pids}" "\${pid}"
        else
          next_pids="\${pid}"
        fi
      else
        if wait "\${pid}"; then
          :
        else
          return \$?
        fi
      fi
    done <<EOF
\${pids}
EOF

    pids="\${next_pids}"
  }

  wait_for_slot() {
    while (( \$(count_pids) >= parallel_jobs )); do
      local before_count=\$(count_pids)
      reap_finished_pids
      if (( \$(count_pids) == before_count )); then
        sleep 1
      fi
    done
  }

  wait_for_all() {
    while (( \$(count_pids) > 0 )); do
      local before_count=\$(count_pids)
      reap_finished_pids
      if (( \$(count_pids) == before_count )); then
        sleep 1
      fi
    done
  }

  if (( parallel_jobs <= 1 )); then
    while IFS=':' read -r prot files; do
      [[ -n "\${prot}" ]] || continue
      run_gwas_qc "\${prot}" "\${files}"
    done < mapping.txt
  else
    pids=""
    while IFS=':' read -r prot files; do
      [[ -n "\${prot}" ]] || continue
      wait_for_slot
      run_gwas_qc "\${prot}" "\${files}" &
      if [[ -n "\${pids}" ]]; then
        printf -v pids '%s\n%s' "\${pids}" "\$!"
      else
        pids="\$!"
      fi
    done < mapping.txt

    wait_for_all
  fi
  """



}
