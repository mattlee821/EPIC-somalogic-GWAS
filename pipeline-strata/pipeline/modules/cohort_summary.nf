process COHORT_SUMMARY {
  label        'low'
  conda        "${params.env_dir}/py_validate"
  publishDir   "${params.outdir}", mode: 'copy'

  input:
  path metrics_list
  path ldsc_summary
  path herit_summary

  output:
  path "all_metrics.tsv", optional: true
  path "metrics_ldsc-herit.tsv", optional: true

  script:
  """
  python ${projectDir}/bin/generate_cohort_summary.py \\
    --metrics_list ${metrics_list} \\
    --ldsc_file ${ldsc_summary} \\
    --herit_file ${herit_summary} \\
    --outdir .
  """

}
