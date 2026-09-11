process PREPARE_PHENOTYPES {
  label        'low'
  conda        "${projectDir}/../env/py_prep"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/002_prepare-phenotypes" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(phenotype_file, stageAs: 'phenotype_input.txt'), path(sample_file), path(covariate_file, stageAs: 'covariate_input.txt'), val(group_column), val(group_values), val(include_proteins), val(covariates)

  output:
  tuple val(study_id), val(group), path("feature_*.pheno"), path("full.pheno"), path("features.manifest"), path("keep_samples.txt")
  path "protein_summary.tsv"

  script:
  """
  # cache-key: prepare_phenotypes explicit strata whitelist (v4)
  python ${projectDir}/bin/prepare_phenotypes.py \\
    --phenotype_file phenotype_input.txt \\
    --sample_file ${sample_file} \\
    --study_id ${study_id} \\
    --group "${group}" \\
    --group_column "${group_column}" \\
    --group_values "${group_values}" \\
    --covariate_file covariate_input.txt \\
    --covariates "${covariates}" \\
    --include_proteins "${include_proteins}" \\
    --outdir .
  """

}
