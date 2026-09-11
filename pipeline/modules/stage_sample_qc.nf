process STAGE_SAMPLE_QC {
  label        'low'
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/004_sample-qc" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(qc_samples)

  output:
  tuple val(study_id), val(group), path("qc_pass_samples.txt")

  script:
  """
  if [[ "\$(readlink -f ${qc_samples})" != "\$(readlink -f qc_pass_samples.txt)" ]]; then
    cp ${qc_samples} qc_pass_samples.txt
  fi
  """
}
