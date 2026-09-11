process PCA {
  label        'medium'
  conda        "${projectDir}/../env/plink2"

  input:
  tuple val(study_id), val(group), path(step1_pgen), path(step1_pvar), path(step1_psam)

  output:
  tuple val(study_id), val(group), path("pca.eigenvec"), path("pca.eigenval")

  script:
  def prefix = step1_pgen.baseName
  """
  bash ${projectDir}/bin/compute_pcs.sh \\
    --pfile ${prefix} \\
    --outdir .
  """
}
