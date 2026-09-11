process LD_PRUNING {
  label        'medium'
  conda        "${params.env_dir}/plink2"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/006_ld-pruning" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), val(pfile_prefix), path(qc_samples), path(qc_variants), val(maf), val(ld_window_kb), val(ld_step), val(ld_r2)

  output:
  tuple val(study_id), val(group), path("step1_input.pgen"), path("step1_input.pvar"), path("step1_input.psam")

  script:
  """
  bash ${projectDir}/bin/ld_pruning.sh \\
    --pfile ${pfile_prefix} \\
    --keep ${qc_samples} \\
    --extract ${qc_variants} \\
    --study ${study_id} \\
    --window ${ld_window_kb} \\
    --step ${ld_step} \\
    --r2 ${ld_r2} \\
    --maf ${maf} \\
    --outdir .
  """
}
