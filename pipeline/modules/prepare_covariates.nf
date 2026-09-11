process PREPARE_COVARIATES {
  label        'low'
  conda        "${projectDir}/../env/py_prep"

  input:
  tuple val(study_id), val(group), path(covariate_file), path(sample_file), val(group_column), val(cases_value), val(include_covariates), val(include_proteins)

  output:
  tuple val(study_id), val(group), path("covariates.cov")

  script:
  """
  # cache-key: prepare_covariates FID=IID + per-study STUDY filter (v3)
  python ${projectDir}/bin/prepare_covariates.py \\
    --covariate_file ${covariate_file} \\
    --sample_file ${sample_file} \\
    --study_id ${study_id} \\
    --group "${group}" \\
    --group_column "${group_column}" \\
    --cases_value "${cases_value}" \\
    --include_covariates "${include_covariates}" \\
    --include_proteins "${include_proteins}" \\
    --outdir .
  """
}
