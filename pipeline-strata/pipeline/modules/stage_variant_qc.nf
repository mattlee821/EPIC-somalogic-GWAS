process STAGE_VARIANT_QC {
  label        'low'
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/005_variant-qc" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(qc_variants)

  output:
  tuple val(study_id), val(group), path("variant_qc_pass.snplist")

  script:
  """
  if [[ "\$(readlink -f ${qc_variants})" != "\$(readlink -f variant_qc_pass.snplist)" ]]; then
    cp ${qc_variants} variant_qc_pass.snplist
  fi
  """
}
