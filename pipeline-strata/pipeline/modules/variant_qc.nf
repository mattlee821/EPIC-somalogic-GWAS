process VARIANT_QC {
  label        'medium'
  conda        "${projectDir}/../env/plink2"

  input:
  tuple val(study_id), val(pfile_prefix), val(maf), val(info_min)

  output:
  tuple val(study_id), path("variant_qc_pass.snplist")

  script:
  """
  bash ${projectDir}/bin/variant_qc.sh \\
    --pfile ${pfile_prefix} \\
    --study ${study_id} \\
    --maf ${maf} \\
    --info_min ${info_min} \\
    --outdir .
  """
}
