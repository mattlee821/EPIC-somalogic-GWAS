process META_GROUP {
  label        'medium'
  conda        "${projectDir}/../env/py_qc"
  publishDir   params.outdir, mode: 'rellink', saveAs: { name ->
    def metalDir = "${study_id}/${protein_id}/meta/010_METAL"
    def qcDir = "${study_id}/${protein_id}/meta/011_QC"
    if (name == 'meta.tsv.gz' || name == 'metrics.tsv') {
      return "${qcDir}/${name}"
    }
    "${metalDir}/${name}"
  }

  input:
  tuple val(study_id), val(protein_id), path(sumstats_files, stageAs: 'meta_inputs??/*')

  output:
  tuple val(study_id), val("meta"), val(protein_id), path("meta.tsv.gz"), path("metrics.tsv"), emit: all_out
  path "meta.tsv.gz", emit: gwas_files
  path "metrics.tsv", emit: metrics_files
  path "meta.tbl.info", optional: true, emit: info_files
  path "meta.log", emit: log_files

  script:
  def configured = (params.metal_binary ?: 'tools/metal/metal').toString()
  def metal_bin = configured.startsWith('/') ? configured : "${projectDir}/${configured}"
  def staged_files = (sumstats_files instanceof java.util.Collection ? sumstats_files : [sumstats_files])
    .collect { it.toString() }
    .join(',')
  """
  export POLARS_MAX_THREADS="${task.cpus}"
  export RAYON_NUM_THREADS="${task.cpus}"
  export OMP_NUM_THREADS="${task.cpus}"
  export OPENBLAS_NUM_THREADS="${task.cpus}"
  export MKL_NUM_THREADS="${task.cpus}"
  export NUMEXPR_MAX_THREADS="${task.cpus}"
  bash ${projectDir}/bin/run_metal.sh \\
    --input_files "${staged_files}" \\
    --protein_id "${protein_id}" \\
    --study "${study_id}" \\
    --group meta \\
    --metal_binary "${metal_bin}" \\
    --outdir .
  """
}
