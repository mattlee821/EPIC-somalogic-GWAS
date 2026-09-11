process REGENIE_STEP0 {
  label        'regenie_step0'
  conda        "${params.env_dir}/regenie"
  publishDir   params.outdir, mode: 'copy', saveAs: { name ->
    if (name == "regenie_step0.master" || name == "step0_split.log" || name ==~ /step0_l0_job\d+\.log/) {
      return "${study_id}/_shared/${group}/006b_regenie-step0/${name}"
    }
    return null
  }

  input:
  tuple val(study_id), val(group), path(step1_pgen), path(step1_pvar), path(step1_psam), path(pheno_file), path(cov_file), val(bsize)

  output:
  tuple val(study_id), val(group), path("regenie_step0.master"), path("regenie_step0_job*.snplist"), path("regenie_step0_job*_l0_*")

  script:
  def pfile_prefix = step1_pgen.baseName
  """
  # REGENIE step 0 must write its master file in the shared work directory
  # because the master records absolute paths to level-0 shard outputs.
  bash ${projectDir}/bin/run_regenie_step0.sh \\
    --step1_pfile ${pfile_prefix} \\
    --pheno_file ${pheno_file} \\
    --cov_file ${cov_file} \\
    --study ${study_id} \\
    --group ${group} \\
    --bsize ${bsize} \\
    --outdir .
  """
}
