process ADD_PCS {
  label        'low'
  conda        "${params.env_dir}/py_validate"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/003_prepare-covariates" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(cov_file, stageAs: 'covariates_in.cov'), path(pcs_file), path(eval_file)

  output:
  tuple val(study_id), val(group), path("covariates.cov")

  script:
  """
  python ${projectDir}/bin/merge_pcs.py \\
    --cov_file ${cov_file} \\
    --pcs_file ${pcs_file} \\
    --eval_file ${eval_file} \\
    --outdir .
  """
}
