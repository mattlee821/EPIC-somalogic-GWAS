#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { VALIDATE_INPUTS } from './modules/validate_inputs.nf'
include { PREPARE_PHENOTYPES } from './modules/prepare_phenotypes.nf'
include { PREPARE_COVARIATES } from './modules/prepare_covariates.nf'
include { VARIANT_QC } from './modules/variant_qc.nf'
include { STAGE_VARIANT_QC } from './modules/stage_variant_qc.nf'
include { LD_PRUNING } from './modules/ld_pruning.nf'
include { PCA } from './modules/pca.nf'
include { ADD_PCS } from './modules/add_pcs.nf'
include { REGENIE_STEP0 } from './modules/regenie_step0.nf'
include { REGENIE_STEP1 } from './modules/regenie_step1.nf'
include { REGENIE_STEP2 } from './modules/regenie_step2.nf'
include { GWAS_QC } from './modules/gwas_qc.nf'
include { META_STUDY } from './modules/meta_study.nf'
include { META_GROUP } from './modules/meta_group.nf'
include { LDSC } from './modules/ldsc.nf'
include { HERITABILITY } from './modules/heritability.nf'
include { COHORT_SUMMARY } from './modules/cohort_summary.nf'

workflow {

    // 0. Resolve and validate the strata definition before anything runs.
    //    A run that silently analyses the wrong samples is far worse than one
    //    that refuses to start, so every failure mode here is fatal.
    if (!(params.strata instanceof Map) || params.strata.isEmpty()) {
        throw new IllegalArgumentException("params.strata is empty. Define at least one stratum in params.yaml.")
    }

    def strata_def = params.strata.collectEntries { name, spec ->
        def column = (spec?.column ?: "").toString().trim()
        def values = (spec?.values ?: []).collect { it.toString().trim() }.findAll { it }
        if (column && !values) {
            throw new IllegalArgumentException(
                "Stratum '${name}' sets column '${column}' but lists no values. " +
                "Membership must be an explicit whitelist -- there is no implicit complement group."
            )
        }
        if (!column && values) {
            throw new IllegalArgumentException(
                "Stratum '${name}' lists values ${values} but no column to test them against."
            )
        }
        [(name.toString()): [column: column, values: values.join(',')]]
    }

    def pca_group = (params.pca_group ?: "").toString().trim()
    if (!strata_def.containsKey(pca_group)) {
        throw new IllegalArgumentException(
            "params.pca_group '${pca_group}' is not a defined stratum. " +
            "Defined: ${strata_def.keySet().sort().join(', ')}"
        )
    }
    if (strata_def[pca_group].column) {
        throw new IllegalArgumentException(
            "params.pca_group '${pca_group}' must be an unfiltered stratum (column: \"\"). " +
            "Ancestry PCs computed inside an exposure-defined subgroup can absorb real " +
            "association signal and make strata non-comparable."
        )
    }

    // analysis_groups holds several names in one CSV field. Nextflow's splitCsv
    // does NOT honour quoted fields unless quote: is passed explicitly -- with
    // it omitted, "all,ever,never" parses as the single value '"all' and every
    // later column shifts. Both splitCsv calls below therefore pass quote, and
    // this parser additionally tolerates a pipe separator and strips any quote
    // that still survives, so a hand-edited samplesheet cannot silently produce
    // the wrong strata.
    def parse_groups = { raw ->
        def text = (raw == null ? "" : raw.toString()).replaceAll('"', '').trim()
        if (!text) return []
        return text.split(/[,|]/).collect { it.trim() }.findAll { it }
    }

    // Every group name the samplesheet requests must key into params.strata.
    // This is checked EAGERLY, by reading the samplesheet here rather than
    // relying on the equivalent check inside the channel .map below: a closure
    // on a channel only throws once that channel is consumed, which is deep
    // into a live run (and never at all under -preview), so a typo'd stratum
    // name could otherwise reach a real submission.
    def samplesheet_groups = file(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .collectMany { row ->
            def parsed = parse_groups(row.analysis_groups)
            parsed ?: params.analysis_groups.collect { it.toString() }
        }
        .unique()

    def unknown_groups = samplesheet_groups.findAll { !strata_def.containsKey(it) }
    if (unknown_groups) {
        throw new IllegalArgumentException(
            "Undefined strata requested: ${unknown_groups.sort().join(', ')}. " +
            "Defined in params.strata: ${strata_def.keySet().sort().join(', ')}"
        )
    }

    // Compact encoding for the validator: "name:column:v1|v2;name2::".
    def strata_spec = strata_def.collect { name, spec ->
        "${name}:${spec.column}:${spec.values.replace(',', '|')}"
    }.join(';')

    // 1. Validation
    // Check if covariate file exists, otherwise fallback to phenotype file for staging
    def cov_file = params.covariate_file && file(params.covariate_file).exists() ? file(params.covariate_file) : file(params.phenotype_file)

    VALIDATE_INPUTS(
        file(params.samplesheet),
        file(params.phenotype_file),
        cov_file,
        params.include_studies,
        params.include_proteins,
        params.covariates,
        params.chromosomes.collect { it.toString() }.join(','),
        strata_spec
    )

    // 2. Load samplesheet and create study channel
    def required_studies = params.include_studies ? params.include_studies.split(',').collect { it.trim() } : []

    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .filter { row -> required_studies.size() == 0 || required_studies.contains(row.study_id) }
        .set { ch_studies }

    def ch_selected_study_count = ch_studies
        .map { row -> row.study_id.toString() }
        .unique()
        .collect()
        .map { study_ids -> study_ids.size() }

    // 3. Prepare covariates and phenotypes, branched by each study's own analysis
    //    groups. The group list comes from the samplesheet 'analysis_groups' column
    //    (set by 004 from per-arm sample counts); it falls back to
    //    params.analysis_groups for samplesheets without that column.
    def ch_study_group = ch_studies
        .flatMap { row ->
            def grps = parse_groups(row.analysis_groups) ?: params.analysis_groups.collect { it.toString() }
            // The PCA stratum is mandatory for every study: PCs are computed on
            // it and shared with all others, so a study that omits it would
            // produce no PCs, join empty at ADD_PCS, and silently emit ZERO
            // results with no error. Force it in rather than trust the
            // samplesheet.
            if (!grps.contains(pca_group)) {
                grps = [pca_group] + grps
            }
            grps.unique().collect { g -> tuple(row, g) }
        }
        .map { row, group ->
            if (!strata_def.containsKey(group)) {
                throw new IllegalArgumentException(
                    "Study ${row.study_id} requests group '${group}', which is not defined in " +
                    "params.strata. Defined: ${strata_def.keySet().sort().join(', ')}"
                )
            }
            tuple(row, group, strata_def[group].column, strata_def[group].values)
        }

    // The stratifying column and its accepted values come from params.strata,
    // not from the samplesheet: the coding of a variable is a property of the
    // data dictionary, identical across studies, so it belongs in one place.
    ch_study_group
        .map { row, group, group_column, group_values ->
            tuple(
                row.study_id,
                group,
                file(params.phenotype_file),
                file(row.sample_file),
                cov_file,
                group_column,
                group_values,
                params.include_proteins,
                params.covariates
            )
        }
        .set { ch_pheno_inputs }

    ch_study_group
        .map { row, group, group_column, group_values ->
            tuple(row.study_id, group, cov_file, file(row.sample_file), group_column, group_values, params.covariates, params.include_proteins)
        }
        .set { ch_cov_inputs }

    PREPARE_PHENOTYPES(ch_pheno_inputs)
    PREPARE_COVARIATES(ch_cov_inputs)

    // 4. Variant QC per study. Sample exclusions (ancestry/relatedness) and HWE
    //    flagging were performed upstream by the imputation pipeline; 002 already
    //    merged the exclusions into the per-study pfile. Variant QC here only
    //    applies the MAF and imputation-quality (R2) filters. (Ancestry PCs are
    //    recomputed on the analysis subset below, not reused from upstream.)
    //    VARIANT_QC takes [sid, pfile, maf, info_min]
    VARIANT_QC(ch_studies.map { row -> tuple(row.study_id, row.pfile, params.maf, params.info_min) })

    // Stage variant QC output into each group directory this study runs
    def ch_sid_group = ch_study_group.map { row, group, group_column, group_values -> tuple(row.study_id, group) }
    STAGE_VARIANT_QC(
        VARIANT_QC.out.combine(ch_sid_group, by: 0)
                     .map { sid, qc_var, grp -> tuple(sid, grp, qc_var) }
    )

    // 4b. Intersect the variant-QC pass set with phenotyped samples for LD pruning
    pheno_samples = PREPARE_PHENOTYPES.out[0].map { sid, grp, feature_files, full, manifest, keep -> tuple(sid, grp, keep) }

    ch_for_ld = ch_studies.map { row -> tuple(row.study_id, row.pfile) }
                          .join(VARIANT_QC.out, by: 0)            // [sid, pfile, qc_variants]
                          .combine(pheno_samples, by: 0)          // [sid, pfile, qc_variants, grp, keep]
                          .map { sid, pfile, qcv, grp, keep ->
                              tuple(sid, grp, pfile, keep, qcv, params.maf, params.ld_window_kb, params.ld_step, params.ld_r2)
                          }

    LD_PRUNING(ch_for_ld)

    // 4c. Recompute ancestry PCs on the analysis subset (proteomics-matched,
    //     deduped, study >= min N). The staged upstream eigenvec/eigenval were
    //     computed on the FULL study, so they are intentionally not used here.
    //     PCs are computed once per study on the params.pca_group LD-pruned set
    //     (the full, unstratified analysis sample) and reused across every
    //     stratum — never computed within a stratum, which could absorb true
    //     association signal and would make strata non-comparable.
    def ch_pca_in = LD_PRUNING.out.filter { sid, grp, pgen, pvar, psam -> grp == pca_group }
    PCA(ch_pca_in)

    def ch_study_pcs = PCA.out.map { sid, grp, pvec, pval -> tuple(sid, pvec, pval) }

    def cov_with_pcs_in = PREPARE_COVARIATES.out
        .combine(ch_study_pcs, by: 0)                            // [sid, grp, cov, eigenvec, eigenval]
        .map { sid, grp, cov, pvec, pval -> tuple(sid, grp, cov, pvec, pval) }

    ADD_PCS(cov_with_pcs_in)

    // 5. Expand to one phenotype file per study, group, and feature
    def pheno_features = PREPARE_PHENOTYPES.out[0]
        .map { sid, grp, feature_files, full, manifest, keep -> tuple(sid, grp, feature_files, manifest) }
        .splitCsv(elem: 3, header: true, sep: '\t')
        .map { row ->
            def protein_id = row[3].protein_id
            def pheno_name = new File(row[3].pheno_file).name
            def feature_files = row[2] instanceof java.util.Collection ? row[2] : [row[2]]
            def pheno_file = feature_files.find { it.name == pheno_name }
            if (pheno_file == null) {
                throw new IllegalStateException("Could not resolve phenotype file '${pheno_name}' for ${row[0]} ${row[1]} ${protein_id}.")
            }
            tuple(row[0], row[1], protein_id, pheno_file)
        }

    // 6. REGENIE Step 0 (shared level-0 work per study and group)
    def regenie_key = { sid, grp -> "${sid}__${grp}" }
    def regenie_full = PREPARE_PHENOTYPES.out[0].map { sid, grp, feature_files, full, manifest, keep -> tuple(regenie_key(sid, grp), full) }
    def regenie_cov = ADD_PCS.out[0].map { sid, grp, cov -> tuple(regenie_key(sid, grp), cov) }
    def regenie_ld = LD_PRUNING.out.map { sid, grp, pgen, pvar, psam -> tuple(regenie_key(sid, grp), sid, grp, pgen, pvar, psam) }

    def regenie_shared = regenie_ld
        .join(regenie_full, by: 0)
        .join(regenie_cov, by: 0)
        .map { key, sid, grp, pgen, pvar, psam, full, cov ->
            tuple(key, sid, grp, pgen, pvar, psam, full, cov)
        }

    step0_inputs = regenie_shared
        .map { key, sid, grp, pgen, pvar, psam, full, cov ->
            tuple(sid, grp, pgen, pvar, psam, full, cov, params.regenie_bsize_step1)
        }

    REGENIE_STEP0(step0_inputs)

    // 7. REGENIE Step 1 (level-1 models per study, group, and feature)
    def pheno_prep = pheno_features.map { sid, grp, protein_id, pheno_file -> tuple(regenie_key(sid, grp), sid, grp, protein_id) }
    def step0_master = REGENIE_STEP0.out.map { sid, grp, master, snplists, l0_files ->
        tuple(regenie_key(sid, grp), master.toString())
    }

    step1_inputs = pheno_prep
        .combine(regenie_shared.map { key, sid, grp, pgen, pvar, psam, full, cov -> tuple(key, pgen, pvar, psam, full, cov) }, by: 0)
        .combine(step0_master, by: 0)
        .map { key, sid, grp, protein_id, pgen, pvar, psam, full, cov, master_file ->
            tuple(sid, grp, protein_id, pgen, pvar, psam, full, cov, master_file, params.regenie_bsize_step1)
        }

    REGENIE_STEP1(step1_inputs)

    // 8. REGENIE Step 2 (one run per study, group, and feature across all chromosomes)
    def chromosome_list = params.chromosomes.collect { it.toString() }.join(',')
    def ch_step1 = REGENIE_STEP1.out.map { sid, grp, protein_id, pred, loco -> tuple("${sid}__${grp}__${protein_id}", pred, loco) }
    def ch_pheno_step2 = pheno_features.map { sid, grp, protein_id, pheno_file -> tuple("${sid}__${grp}__${protein_id}", sid, grp, protein_id, pheno_file) }
    def ch_cov_step2 = ADD_PCS.out[0].map { sid, grp, cov -> tuple("${sid}_${grp}", cov) }
    def ch_genotype = ch_studies.map { r -> tuple(r.study_id, r.pfile, r.sample_file) }

    step2_setup = ch_pheno_step2
        .join(ch_step1, by: 0)
        .map { key, sid, grp, protein_id, pheno_file, pred, loco -> tuple("${sid}_${grp}", sid, grp, protein_id, pheno_file, pred, loco) }
        .combine(ch_cov_step2, by: 0)
        .map { key, sid, grp, protein_id, pheno_file, pred, loco, cov -> tuple(sid, grp, protein_id, pheno_file, pred, loco, cov) }
        .combine(ch_genotype, by: 0)
        .map { sid, grp, protein_id, pheno_file, pred, loco, cov, pfile, samp ->
            tuple(sid, grp, protein_id, pfile, samp, chromosome_list, pheno_file, cov, pred, loco, params.regenie_bsize_step2)
        }

    REGENIE_STEP2(step2_setup)

    // 9. Group REGENIE step 2 outputs by protein
    ch_grouped_regenie = REGENIE_STEP2.out
        .map { sid, grp, protein_id, files ->
            def result_files = files instanceof java.util.Collection ? files : [files]
            tuple(sid, grp, protein_id, result_files)
        }

    // 10. Consolidated post-GWAS QC
    def gwas_qc_batch_size = (params.gwas_qc_batch_size ?: 50) as int
    if (gwas_qc_batch_size < 1) {
        throw new IllegalArgumentException("params.gwas_qc_batch_size must be >= 1")
    }

    // Batch deterministically: groupTuple(size:) chunks on SLURM *arrival* order,
    // so both batch membership and within-batch ordering vary run-to-run. That
    // changes the GWAS_QC task hash (prot_ids/files/mapping) every run, so it
    // never caches even though all STEP2 inputs are cached. Instead collect the
    // full per-(study,group) group, sort by protein_id, then slice into fixed
    // chunks -> identical batches every run -> cacheable on -resume.
    ch_gwas_qc_in = ch_grouped_regenie
        .map { sid, grp, prot, fs ->
            def files = (fs instanceof java.util.Collection ? fs : [fs])
            tuple("${sid}__${grp}", sid, grp, prot.toString(), files)
        }
        .groupTuple(by: 0)
        .flatMap { key, sids, grps, prots, file_sets ->
            def rows = (0..<prots.size())
                .collect { i -> [prots[i].toString(), file_sets[i], sids[i].toString(), grps[i].toString()] }
                .sort { a, b -> a[0] <=> b[0] }
            rows.collate(gwas_qc_batch_size).collect { batch ->
                def prot_ids = batch.collect { it[0] }
                def files = batch.collectMany { it[1] instanceof java.util.Collection ? it[1] : [it[1]] }
                def mapping = batch.collect { row ->
                    def prot_files = row[1] instanceof java.util.Collection ? row[1] : [row[1]]
                    "${row[0]}:${prot_files.collect { it.name }.join(',')}"
                }.join('\n')
                tuple(batch[0][2], batch[0][3], prot_ids, files, mapping)
            }
        }

    GWAS_QC(ch_gwas_qc_in)

    // Flatten batched GWAS_QC outputs for downstream use
    ch_gwas_qc_results = GWAS_QC.out.all_out
        .flatMap { sid, grp, prots, gwas_files, metrics_files ->
            def prot_ids = prots instanceof java.util.Collection ? prots.collect { it.toString() } : [prots.toString()]
            def gwas_list = gwas_files instanceof java.util.Collection ? gwas_files : [gwas_files]
            def metrics_list = metrics_files instanceof java.util.Collection ? metrics_files : [metrics_files]

            def protein_from_output = { f ->
                def parent = f.getParent()?.getFileName()?.toString()
                if (parent) {
                    return parent
                }
                def parts = f.toString().tokenize('/\\')
                return parts.size() >= 2 ? parts[-2] : null
            }

            def gwas_by_prot = gwas_list.collectEntries { f ->
                [(protein_from_output(f)): f]
            }
            def metrics_by_prot = metrics_list.collectEntries { f ->
                [(protein_from_output(f)): f]
            }

            prot_ids.collect { prot_id ->
                def gwas = gwas_by_prot[prot_id]
                def metrics = metrics_by_prot[prot_id]
                if (gwas == null || metrics == null) {
                    throw new IllegalStateException(
                        "Could not resolve GWAS_QC outputs for ${sid} ${grp} ${prot_id}. " +
                        "GWAS keys=${gwas_by_prot.keySet().findAll { it != null }.sort().join(',')} " +
                        "Metrics keys=${metrics_by_prot.keySet().findAll { it != null }.sort().join(',')}"
                    )
                }
                tuple(sid, grp, prot_id, gwas, metrics)
            }
        }

    // 8b. Within-study METAL meta across a PAIR of strata (optional).
    //     params.meta_group names the two strata, e.g. "ever,never". For an
    //     exposure split this is usually NOT what you want: pooling ever+never
    //     just reconstructs the unstratified analysis. Compare per-stratum
    //     effect sizes instead. Kept because it is occasionally useful for
    //     sensitivity checks, and disabled by default.
    def meta_group_raw = (params.meta_group ?: "").toString().trim()
    def meta_group_names = (meta_group_raw.toLowerCase() in ['', 'false', '0', 'no']) ?
        [] : meta_group_raw.split(',').collect { it.trim() }.findAll { it }
    if (meta_group_names && meta_group_names.size() != 2) {
        throw new IllegalArgumentException(
            "params.meta_group must name exactly two strata (e.g. \"ever,never\") or be empty; got: ${meta_group_names}"
        )
    }
    meta_group_names.each { name ->
        if (!strata_def.containsKey(name)) {
            throw new IllegalArgumentException(
                "params.meta_group references '${name}', which is not a defined stratum. " +
                "Defined: ${strata_def.keySet().sort().join(', ')}"
            )
        }
    }
    def meta_group_enabled = meta_group_names.size() == 2
    def ch_meta_group_qc_out = Channel.empty()
    def ch_meta_group_for_meta = Channel.empty()

    if (meta_group_enabled) {
        def ch_meta_group_in = ch_gwas_qc_results
            .filter { sid, grp, prot, gwas, metrics -> meta_group_names.contains(grp.toString()) }
            .map { sid, grp, prot, gwas, metrics -> tuple("${sid}__${prot}", sid, prot, gwas) }
            .groupTuple(by: 0, size: 2, remainder: true)
            .filter { key, sids, prots, gwases -> gwases.size() == 2 }
            .map { key, sids, prots, gwases ->
                def files = gwases instanceof java.util.Collection ? gwases : [gwases]
                tuple(sids[0].toString(), prots[0].toString(), files)
            }

        META_GROUP(ch_meta_group_in)

        ch_meta_group_qc_out = META_GROUP.out.all_out
        ch_meta_group_for_meta = META_GROUP.out.all_out
            .map { sid, grp, prot, gwas, met -> tuple(sid, grp, prot, gwas) }
    }

    // 9. Meta-analysis across studies, run independently within each stratum
    //    named by params.meta_study, so each stratum gets its own pooled sumstats.
    def meta_study_raw = params.meta_study
    def meta_study_enabled = meta_study_raw && !(meta_study_raw.toString().trim().toLowerCase() in ['false', '0', 'no'])

    // Per-protein LDSC / heritability rows, merged to one file each (or null when
    // disabled). COHORT_SUMMARY joins them into metrics_ldsc-herit.tsv.
    def ch_ldsc_summary = null
    def ch_herit_summary = null

    if (meta_study_enabled) {
        def meta_study_list = meta_study_raw.toString().split(',').collect { it.trim().toLowerCase() }.findAll { it }
        if (meta_study_list == ["true"] || meta_study_list == ["1"]) meta_study_list = ["meta"]
        def meta_study_set = meta_study_list as Set

        ch_gwas_qc_results.map { sid, grp, prot, gwas, met -> tuple(sid, grp, prot, gwas) }
            .mix(ch_meta_group_for_meta)
            .filter { sid, grp, prot, gwas -> meta_study_set.contains(grp.toString().toLowerCase()) }
            .combine(ch_selected_study_count)
            .map { sid, grp, prot, gwas, study_count ->
                tuple(groupKey("${grp}__${prot}", study_count as int), sid, grp, prot, gwas)
            }
            .groupTuple(remainder: true)
            .filter { key, sids, grps, prots, gwases -> gwases.size() > 1 }
            .map { key, sids, grps, prots, gwases ->
                def files = gwases instanceof java.util.Collection ? gwases : [gwases]
                tuple(grps[0].toString(), prots[0].toString(), files)
            }
            .set { ch_meta_study_in }

        META_STUDY(ch_meta_study_in)

        // 9b. LDSC intercept + SNP-h2 on each cross-study meta sumstats.
        //     Diagnostic: intercept ~1 => the meta lambda inflation is polygenic,
        //     not confounding/stratification. Runs per protein; rows are merged
        //     and handed to COHORT_SUMMARY for the combined metrics_ldsc-herit.tsv.
        def ldsc_enabled = params.ldsc && !(params.ldsc.toString().trim().toLowerCase() in ['false', '0', 'no'])
        if (ldsc_enabled) {
            LDSC(
                META_STUDY.out.all_out.map { study, grp, prot, gwas, met -> tuple(grp, prot, gwas) }
            )
            ch_ldsc_summary = LDSC.out.row.collectFile(name: 'ldsc_rows.tsv', keepHeader: true, skip: 1, sort: true)
        }

        // 9c. Heritability decomposition (UKB-PPP style): per-protein total SNP-h2
        //     = pQTL component (cis/trans variance from clumped index pQTLs)
        //     + polygenic component (LDSC --h2 on the region-excluded genome).
        def herit_enabled = params.herit && !(params.herit.toString().trim().toLowerCase() in ['false', '0', 'no'])
        if (herit_enabled) {
            HERITABILITY(
                META_STUDY.out.all_out.map { study, grp, prot, gwas, met -> tuple(grp, prot, gwas) }
            )
            ch_herit_summary = HERITABILITY.out.row.collectFile(name: 'herit_rows.tsv', keepHeader: true, skip: 1, sort: true)
        }
    }

    // 10. Cohort summary
    ch_metrics_all = ch_gwas_qc_results.map { sid, grp, prot, gwas, met -> met.toString() }

    if (meta_group_enabled) {
        ch_metrics_all = ch_metrics_all.mix(ch_meta_group_qc_out.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    if (meta_study_enabled) {
        ch_metrics_all = ch_metrics_all.mix(META_STUDY.out.all_out.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    def ch_metrics_list = ch_metrics_all.collectFile(name: 'metrics_list.txt', newLine: true)

    // LDSC/herit are optional: pass an empty placeholder when a side is disabled
    // so COHORT_SUMMARY still runs and writes all_metrics.tsv.
    def no_file = file("${projectDir}/assets/NO_FILE")
    COHORT_SUMMARY(
        ch_metrics_list,
        ch_ldsc_summary  ?: Channel.value(no_file),
        ch_herit_summary ?: Channel.value(no_file)
    )
}
