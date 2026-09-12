process META_STUDY {
  label        'medium'
  conda        "${params.env_dir}/py_qc"
  publishDir   params.outdir, mode: 'copy', overwrite: true, saveAs: { name ->
    // The stratum MUST be in the path. Several strata are meta-analysed in one
    // run, and with overwrite:true a group-less path makes them silently
    // clobber each other -- the last stratum to finish would win, for every
    // protein. Layout matches what src/006_staging.sh already expects.
    def metalDir = "meta/${protein_id}/${group}/012_METAL"
    def qcDir = "meta/${protein_id}/${group}/013_QC"
    if (name == 'meta.tsv.gz' || name == 'metrics.tsv') {
      return "${qcDir}/${name}"
    }
    "${metalDir}/${name}"
  }

  input:
  tuple val(group), val(protein_id), path(sumstats_files, stageAs: 'meta_inputs??/*')

  output:
  tuple val("meta_study"), val(group), val(protein_id), path("meta.tsv.gz"), path("metrics.tsv"), emit: all_out
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
    --study meta_study \\
    --group "${group}" \\
    --metal_binary "${metal_bin}" \\
    --outdir .
  """
}
