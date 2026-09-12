process LDSC {
  label        'low'
  conda        "${params.env_dir}/ldsc"
  publishDir   params.outdir, mode: 'copy', overwrite: true, saveAs: { name ->
    // Per-protein metrics + logs alongside the meta outputs.
    "meta/${protein_id}/${group}/014_metrics/${file(name).name}"
  }

  input:
  tuple val(group), val(protein_id), path(meta_tsv)

  output:
  path "metrics_ldsc.tsv", emit: row
  path "ldsc.run.log", emit: log
  path "h2.log", optional: true, emit: h2log

  script:
  """
  bash ${projectDir}/bin/run_ldsc.sh \\
    --meta "${meta_tsv}" \\
    --protein_id "${protein_id}" \\
    --ldsc_home "${params.ldsc_home}" \\
    --ref_dir "${params.ldsc_ref_dir}" \\
    --hm3 "${params.ldsc_hm3}" \\
    --outdir . \\
    --out_row metrics_ldsc.tsv
  """
}
