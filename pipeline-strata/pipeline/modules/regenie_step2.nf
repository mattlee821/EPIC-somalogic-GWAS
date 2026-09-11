process REGENIE_STEP2 {
  label        'regenie_step2'
  conda        "${params.env_dir}/regenie"
  publishDir   params.outdir, mode: params.publish_mode ?: 'link', saveAs: { name ->
    // Single all-chromosome output per protein (full_<protein>.regenie.gz).
    "${study_id}/${protein_id}/${group}/008_regenie-step2/${file(name).name}"
  }

  input:
  tuple val(study_id), val(group), val(protein_id), val(pfile), val(sample_file), val(chromosomes), path(pheno_file), path(cov_file), path(pred_list), path(loco_files), val(bsize)

  output:
  tuple val(study_id), val(group), val(protein_id), path("*.regenie.gz")

  script:
  """
  # cache-key: step2 single-pass all-chromosome run (v2)
  bash ${projectDir}/bin/run_regenie_step2.sh \\
    --pfile ${pfile} \\
    --sample_file ${sample_file} \\
    --pheno_file ${pheno_file} \\
    --cov_file ${cov_file} \\
    --pred_list ${pred_list} \\
    --study ${study_id} \\
    --group ${group} \\
    --chromosomes "${chromosomes}" \\
    --bsize ${bsize} \\
    --outdir .
  """
}
