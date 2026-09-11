process REGENIE_STEP1 {
  label        'regenie_step1'
  conda        "${projectDir}/../env/regenie"
  publishDir   { "${params.outdir}/${study_id}/${protein_id}/${group}/007_regenie-step1" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), val(protein_id), path(step1_pgen), path(step1_pvar), path(step1_psam), path(pheno_file), path(cov_file), val(master_file), val(bsize)

  output:
  tuple val(study_id), val(group), val(protein_id), path("pred.list"), path("loco_*.loco.gz")

  script:
  def pfile_prefix = step1_pgen.baseName
  """
  # Refresh task hash after pred.list localization fix so -resume reruns Step 1.
  bash ${projectDir}/bin/run_regenie_step1.sh \\
    --step1_pfile ${pfile_prefix} \\
    --pheno_file ${pheno_file} \\
    --cov_file ${cov_file} \\
    --master_file "${master_file}" \\
    --phenotype_id "${protein_id}" \\
    --study ${study_id} \\
    --group ${group} \\
    --bsize ${bsize} \\
    --outdir .
  """

}
