process PREPARE_COVARIATES {
  label        'low'
  conda        "${params.env_dir}/py_prep"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/002_prepare-covariates" }, mode: 'copy', pattern: 'covariates_used.txt'

  input:
  tuple val(study_id), val(group), path(covariate_file), path(sample_file), val(group_column), val(group_values), val(include_covariates), val(include_proteins)

  output:
  tuple val(study_id), val(group), path("covariates.cov"), emit: cov
  path "covariates_used.txt", emit: used

  script:
  """
  # cache-key: prepare_covariates explicit strata whitelist + psam covariate fallback (v4)
  python ${projectDir}/bin/prepare_covariates.py \\
    --covariate_file ${covariate_file} \\
    --sample_file ${sample_file} \\
    --study_id ${study_id} \\
    --group "${group}" \\
    --group_column "${group_column}" \\
    --group_values "${group_values}" \\
    --include_covariates "${include_covariates}" \\
    --include_proteins "${include_proteins}" \\
    --outdir .
  """
}
