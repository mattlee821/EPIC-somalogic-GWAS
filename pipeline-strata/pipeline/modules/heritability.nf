process HERITABILITY {
  label        'low'
  conda        "${projectDir}/../env/herit"
  publishDir   params.outdir, mode: 'copy', overwrite: true, saveAs: { name ->
    "meta/${protein_id}/${group}/014_metrics/${file(name).name}"
  }

  input:
  tuple val(group), val(protein_id), path(meta_tsv)

  output:
  path "metrics_herit.tsv", emit: row
  path "herit.run.log", emit: log

  script:
  def ld_ref = params.ld_ref_pfile ? "--ld_ref \"${params.ld_ref_pfile}\"" : ""
  """
  bash ${projectDir}/bin/run_herit.sh \\
    --meta "${meta_tsv}" \\
    --protein_id "${protein_id}" \\
    --protein_map "${params.protein_gene_map}" \\
    ${ld_ref} \\
    --ldsc_home "${params.ldsc_home}" \\
    --ref_dir "${params.ldsc_ref_dir}" \\
    --hm3 "${params.ldsc_hm3}" \\
    --cis_window_kb "${params.cis_window_kb}" \\
    --clump_p "${params.herit_p}" \\
    --clump_r2 "${params.herit_clump_r2}" \\
    --hla "${params.hla_region}" \\
    --outdir . \\
    --out_row metrics_herit.tsv
  """
}
