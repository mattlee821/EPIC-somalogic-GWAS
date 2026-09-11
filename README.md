# EPIC Somalogic GWAS Pipeline

**Version:** May 2026

This repository contains a Nextflow DSL2 pipeline for protein GWAS across EPIC
studies. The current pipeline runs study/group-specific sample and variant QC,
phenotype/covariate preparation, PCA, REGENIE Step 0/1/2, post-GWAS QC, optional
within-study and across-study METAL meta-analysis, and a final cohort metrics
summary.

The finalized workflow lives in `pipeline/`. The helper scripts used to prepare
local inputs and submit the pipeline live in `src/`.

## Quick Start

### 0. Create `.env`

Start by creating a project `.env` file from the template and editing the paths
for the current system:

```bash
cp .env.example .env
```

At minimum, check these values before running any scripts:

- `GWAS_PROJECT_ROOT`
- `GWAS_ANALYSIS_DIR`
- `GWAS_GENETICS_DIR`
- `GWAS_PHENO_FILE`
- `GWAS_SAMPLESHEET`
- `GWAS_EXTERNAL_SOMA_RDS`
- `GWAS_DEPOT_GENETICS_DIR`

The `src/00*` helper scripts source `.env` and use those paths to prepare inputs,
run Nextflow, and stage final transfer files.

### 1. Set Up the Runtime and Tools

Run `001_source.sh` to create the `nf_master` environment and install the bundled
METAL binary used by the pipeline:

```bash
bash src/001_source.sh
```

### 2. Prepare Genetics, Phenotypes, and the Samplesheet

Run the input preparation scripts in order:

```bash
bash src/002_data-genetics.sh
Rscript src/003_data-phenotype.R
Rscript src/004_samplesheet.R
```

These create or update the genetics layout, phenotype files, trait lists, and
`pipeline/samplesheet.csv`.

### 3. Run the GWAS Pipeline

Run a small test through the SLURM wrapper first:

```bash
sbatch src/005_run.sh \
  --profile slurm \
  --out analysis/test-1000 \
  --studies "Neuro_01,Brea_01_Erneg" \
  --proteins "data/phenotype/traits-1000.txt" \
  --pheno "data/phenotype/phenofile.txt" \
  --covariates "Sex,Age_Recr" \
  --meta-group true \
  --meta-study "meta"
```

Then run all numeric protein traits by passing an empty protein selector:

```bash
sbatch src/005_run.sh \
  --out analysis/full \
  --studies "Neuro_01,Brea_01_Erneg" \
  --proteins "" \
  --pheno "data/phenotype/phenofile.txt" \
  --covariates "Sex,Age_Recr" \
  --meta-group true \
  --meta-study "meta"
```

`src/005_run.sh` always runs Nextflow with `-resume`. Edit the `#SBATCH` header
at the top of that file to change the submit job name, logs, walltime, memory,
CPUs, or partition.

### 4. Stage Transfer-Ready Results

After the pipeline and across-study meta-analysis have completed, run
`006_staging.sh` to collect final `meta.tsv.gz` files into a transfer directory:

```bash
bash src/006_staging.sh \
  --out analysis/full \
  --transfer transfer
```

Use `--dry-run` first to inspect what will be staged:

```bash
bash src/006_staging.sh --out analysis/full --transfer transfer --dry-run
```

The staging script writes one transfer file per feature, named
`<feature>.tsv.gz`, plus a staging manifest in the transfer directory. It can
also be submitted through SLURM:

```bash
sbatch src/006_staging.sh --out analysis/full --transfer transfer
```

### Complete `src/00*` Order

For a fresh run, use the numbered helper scripts in this order:

1. Create and edit `.env` from `.env.example`.
2. `bash src/001_source.sh`
3. `bash src/002_data-genetics.sh`
4. `Rscript src/003_data-phenotype.R`
5. `Rscript src/004_samplesheet.R`
6. `sbatch src/005_run.sh ...`
7. `bash src/006_staging.sh --out <analysis-dir> --transfer <transfer-dir>`

The final staging step should be run only after the desired across-study
meta-analysis outputs are present under the analysis directory.

## Inputs

### `pipeline/samplesheet.csv`

Required columns:

- `study_id`
- `pfile`: prefix to one all-chromosome PLINK2 dataset. The pipeline requires
  `<prefix>.pgen`, `<prefix>.pvar`, and `<prefix>.psam`.
- `sample_file`: sample sheet for the study. This must be exactly
  `<pfile>.psam`.
- `group_column`: sample or covariate column used for split analyses.
- `cases_value`: value in `group_column` interpreted as cases.

Example `samplesheet.csv` contents:

| study_id | pfile | sample_file | group_column | cases_value |
|---|---|---|---|---|
| `Study_01` | `../my/file/path/genetics/Study_01/study_01_genotypes` | `../my/file/path/genetics/Study_01/study_01_genotypes.psam` | `PHENO` | `2` |
| `Study_02` | `../my/file/path/genetics/Study_02/study_02_genotypes` | `../my/file/path/genetics/Study_02/study_02_genotypes.psam` | `PHENO` | `2` |

The default analysis groups are defined in `pipeline/params.yaml`:

```yaml
analysis_groups: ["cases", "controls"]
```

The code also supports a `combined` group if it is added to `analysis_groups`.
For `combined`, the group column is retained as a covariate when available.

### Phenotype and Covariate Files

Set these in `pipeline/params.yaml` or override them at runtime:

- `phenotype_file`
- `covariate_file`

If the configured covariate file does not exist, the pipeline falls back to the
phenotype file for covariate staging. Both phenotype and covariate inputs must
contain `IID`. `FID` is added as `0` during preparation if needed.

### Protein and Study Selection

Runtime selectors:

- `--include_studies`: comma-separated study IDs, or empty for all studies in
  the samplesheet.
- `--include_proteins`: comma-separated protein IDs, a path to a one-ID-per-line
  file, or empty for all numeric phenotype columns.
- `--covariates`: comma-separated covariate names. If empty, numeric covariates
  are selected automatically after excluding IDs, metadata, and requested protein
  columns.

The SLURM wrapper defaults to a small test list:

```bash
INCLUDE_PROTEINS="data/phenotype/traits-10.txt"
```

Pass `--proteins ""` to the wrapper when you want all numeric traits.

## Running the Pipeline

The recommended entry point is:

```bash
sbatch src/005_run.sh [options]
```

Supported wrapper options:

- `--profile`: Nextflow profile. Default: `slurm`.
- `--out`, `-o`, `--outdir`: output directory. Default: `<project>/analysis`.
- `--studies`: comma-separated study IDs.
- `--proteins`: comma-separated protein IDs, a list file, or `""` for all
  numeric traits.
- `--pheno`: phenotype file path.
- `--covar`: covariate file path.
- `--covariates`: comma-separated covariate names.
- `--meta-group`: `true` or `false`. Runs within-study cases+controls
  meta-analysis when enabled.
- `--meta-study`: `false`, `true`, or a comma-separated set of analysis groups
  to meta-analyze across studies. `true` maps to `meta`, meaning the
  within-study meta outputs are meta-analyzed across studies. Explicit values can
  include `cases`, `controls`, `combined`, and `meta` if those outputs exist.

The wrapper sources `.env`, sets `GWAS_PROJECT_ROOT`, activates the `nf_master`
conda environment, and runs:

```bash
cd pipeline
nextflow run main.nf -profile <profile> -resume -params-file params.yaml ...
```

## Workflow in `src/`

Primary helper scripts:

- `src/001_source.sh`: create the `nf_master` environment and install the bundled
  METAL binary into `pipeline/tools/metal/metal`.
- `src/002_data-genetics.sh`: project-specific helper to copy documented genetic
  inputs into the local data layout.
- `src/003_data-phenotype.R`: build `data/phenotype/phenofile.txt` and trait list
  files such as `traits-10.txt`, `traits-100.txt`, and `traits-1000.txt`.
- `src/004_samplesheet.R`: generate `pipeline/samplesheet.csv` from the phenotype
  file and genetics directory.
- `src/005_run.sh`: SLURM submission wrapper for the Nextflow pipeline.
- `src/006_staging.sh`: final staging step that collects across-study
  meta-analysis outputs into a transfer/upload directory.

## Pipeline Parameters

Important defaults in `pipeline/params.yaml`:

| Parameter | Default | Meaning |
|---|---:|---|
| `analysis_groups` | `["cases", "controls"]` | Study groups to run |
| `gwas_qc_batch_size` | `10` | Proteins handled per `GWAS_QC` task within each study/group |
| `maf` | `0.01` | Variant MAF threshold |
| `hwe` | `1e-15` | Variant HWE threshold |
| `geno` | `0.05` | Variant missingness threshold |
| `mind` | `0.05` | Sample missingness threshold |
| `info_score` | `0.3` | Minimum INFO score in GWAS QC |
| `king_cutoff` | `0.0884` | KING relatedness cutoff |
| `ld_window_kb` | `1000` | LD-pruning window |
| `ld_step` | `1` | LD-pruning step |
| `ld_r2` | `0.5` | LD-pruning r2 threshold |
| `regenie_bsize_step1` | `1000` | REGENIE Step 0/1 block size |
| `regenie_bsize_step2` | `400` | REGENIE Step 2 block size |
| `chromosomes` | `1..22` | Chromosomes processed in Step 2 |

## Pipeline Steps

### 1. Input Validation

**Module:** `pipeline/modules/validate_inputs.nf`  
**Script:** `pipeline/bin/validate_inputs.py`

Validates samplesheet files, phenotype/covariate files, requested studies,
requested proteins, requested covariates, and group columns. Writes:

```text
<outdir>/_shared/001_validation/validation_report.txt
```

### 2. Phenotype and Covariate Preparation

**Modules:** `prepare_phenotypes.nf`, `prepare_covariates.nf`  
**Scripts:** `prepare_phenotypes.py`, `prepare_covariates.py`

For each selected study and analysis group, the pipeline:

- intersects phenotype/covariate samples with the study sample file;
- filters cases or controls using `group_column` and `cases_value`;
- selects requested proteins, or all numeric non-covariate traits;
- writes one `feature_*.pheno` file per protein plus `full.pheno`;
- writes `features.manifest`, `protein_summary.tsv`, and `keep_samples.txt`;
- prepares a REGENIE covariate file;
- later merges 10 PCs into `covariates.cov`.

### 3. Sample QC

**Module:** `sample_qc.nf`  
**Script:** `sample_qc.sh`

PLINK2 sample filters:

- `--mind 0.05`
- heterozygosity outliers using F coefficient mean +/- 3 SD, computed after
  temporary LD pruning with `--indep-pairwise 200 50 0.25`
- `--king-cutoff 0.0884`

Output staged per group:

```text
<outdir>/<study>/_shared/<group>/004_sample-qc/qc_pass_samples.txt
```

### 4. Variant QC

**Module:** `variant_qc.nf`  
**Script:** `variant_qc.sh`

PLINK2 variant filters:

- `--maf 0.01`
- `--hwe 1e-15`
- `--geno 0.05`

Output staged per group:

```text
<outdir>/<study>/_shared/<group>/005_variant-qc/variant_qc_pass.snplist
```

### 5. LD Pruning and PCA

**Modules:** `ld_pruning.nf`, `pca.nf`, `add_pcs.nf`  
**Scripts:** `ld_pruning.sh`, `compute_pcs.sh`, `merge_pcs.py`

LD pruning uses:

```text
--indep-pairwise 1000kb 1 0.5
```

The resulting `step1_input.*` PLINK2 pfile set is used for REGENIE Step 0/1 and
for PCA. PCA uses `plink2 --pca 10`; the PCs are merged into the group covariate
file.

### 6. REGENIE Step 0

**Module:** `regenie_step0.nf`  
**Script:** `run_regenie_step0.sh`

Step 0 is shared per study and group. It runs REGENIE Step 1 split-level-0 work
with:

- `--step 1`
- `--bsize <regenie_bsize_step1>`
- `--split-l0`
- `--run-l0`
- `--gz`
- `--force-step1`

Published files include the level-0 master file and logs under:

```text
<outdir>/<study>/_shared/<group>/006b_regenie-step0/
```

### 7. REGENIE Step 1

**Module:** `regenie_step1.nf`  
**Script:** `run_regenie_step1.sh`

Step 1 runs one task per study, group, and protein using the Step 0 master file:

- `--step 1`
- `--run-l1 <master>`
- `--l1-phenoList <protein_id>`
- `--keep-l0`
- `--bsize <regenie_bsize_step1>`
- `--gz`
- `--force-step1`

Outputs:

```text
<outdir>/<study>/<protein>/<group>/007_regenie-step1/pred.list
<outdir>/<study>/<protein>/<group>/007_regenie-step1/loco_*.loco.gz
```

### 8. REGENIE Step 2

**Module:** `regenie_step2.nf`  
**Script:** `run_regenie_step2.sh`

Step 2 runs one task per study, group, and protein, then loops over the configured
chromosome list inside that task.

Key settings:

- `--step 2`
- `--pgen <pfile>`
- `--qt`
- `--gz`
- `--bsize <regenie_bsize_step2>`
- `--chr <chromosome>`

Step 2 reads the same all-chromosome PLINK2 prefix for every configured
chromosome and uses `--chr <chromosome>` to write chromosome-specific outputs.

Outputs:

```text
<outdir>/<study>/<protein>/<group>/008_regenie-step2/chr<chromosome>.regenie.gz
```

### 9. Post-GWAS QC

**Module:** `gwas_qc.nf`  
**Script:** `gwas_qc.py`

`GWAS_QC` batches proteins by `gwas_qc_batch_size` within each study/group. For
each protein it:

- concatenates REGENIE chromosome outputs;
- computes `P = 10^(-LOG10P)`;
- renames REGENIE columns to the standard output names;
- filters to `INFO >= info_score`;
- filters to `N >= 0.5 * max(N)` for that protein/study/group;
- computes genomic inflation metrics before and after filtering.

GWAS output columns:

```text
CHR, POS, ID, EA, OA, EAF, BETA, SE, P, N, INFO
```

Outputs:

```text
<outdir>/<study>/<protein>/<group>/009_QC/gwas.tsv.gz
<outdir>/<study>/<protein>/<group>/009_QC/metrics.tsv
```

### 10. Meta-Analysis

**Modules:** `meta_group.nf`, `meta_study.nf`  
**Scripts:** `run_metal.sh`, `metal_qc.py`

METAL settings:

- `TRACKPOSITIONS ON`
- `SCHEME STDERR`
- `AVERAGEFREQ ON`
- `MINMAXFREQ ON`
- `CUSTOMVARIABLE TotalSampleSize`
- `LABEL TotalSampleSize as N`
- `MARKER ID`
- `ALLELE EA OA`
- `EFFECT BETA`
- `STDERR SE`
- `PVALUE P`
- `WEIGHT N`
- `FREQ EAF`
- `CHROMOSOMELABEL CHR`
- `POSITIONLABEL POS`

Within-study meta-analysis is enabled with `--meta-group true`. It combines
cases and controls within each study/protein and writes:

```text
<outdir>/<study>/<protein>/meta/010_METAL/meta.log
<outdir>/<study>/<protein>/meta/010_METAL/meta.tbl.info
<outdir>/<study>/<protein>/meta/011_QC/meta.tsv.gz
<outdir>/<study>/<protein>/meta/011_QC/metrics.tsv
```

Across-study meta-analysis is enabled with `--meta-study`. `--meta-study true`
means `meta`, which meta-analyzes the within-study meta outputs across studies.
Explicit values such as `cases,controls,meta` can be used when those group
outputs exist.

Across-study outputs:

```text
<outdir>/meta/<protein>/<group>/012_METAL/meta.log
<outdir>/meta/<protein>/<group>/012_METAL/meta.tbl.info
<outdir>/meta/<protein>/<group>/013_QC/meta.tsv.gz
<outdir>/meta/<protein>/<group>/013_QC/metrics.tsv
```

Meta QC keeps the standard columns:

```text
CHR, POS, ID, EA, OA, EAF, BETA, SE, P, N
```

It filters to `N >= 0.5 * max(N)` and writes lambda metrics.

### 11. Cohort Summary

**Module:** `cohort_summary.nf`  
**Script:** `generate_cohort_summary.py`

Collects all `metrics.tsv` files from GWAS QC and enabled meta-analysis stages,
drops internal flag columns if present, and writes:

```text
<outdir>/all_metrics.tsv
```

## Output Structure

Representative output layout:

```text
<outdir>/
+-- _shared/
|   +-- 001_validation/
+-- <study>/
|   +-- _shared/
|   |   +-- <group>/
|   |       +-- 002_prepare-phenotypes/
|   |       +-- 003_prepare-covariates/
|   |       +-- 004_sample-qc/
|   |       +-- 005_variant-qc/
|   |       +-- 006_ld-pruning/
|   |       +-- 006b_regenie-step0/
|   +-- <protein>/
|       +-- <group>/
|       |   +-- 007_regenie-step1/
|       |   +-- 008_regenie-step2/
|       |   +-- 009_QC/
|       +-- meta/
|           +-- 010_METAL/
|           +-- 011_QC/
+-- meta/
|   +-- <protein>/
|       +-- <group>/
|           +-- 012_METAL/
|           +-- 013_QC/
+-- all_metrics.tsv
+-- pipeline_info/
    +-- execution_report.html
    +-- execution_timeline.html
    +-- execution_trace.txt
```

## Runtime and Retry Settings

The SLURM profile in `pipeline/nextflow.config` uses a queue size of `1000` and
retries transient signal-like exits (`130`, `137`, `139`, `140`, `143`) up to two
times for most tasks. `REGENIE_STEP0` and `REGENIE_STEP1` terminate immediately
on failure.

Important process time limits:

- `REGENIE_STEP2`: `2.h`
- `GWAS_QC`: `2.h`
- `META_GROUP`: `2.h`
- `META_STUDY`: `2.h`

The explicit `2.h` limits for `GWAS_QC`, `META_GROUP`, and `META_STUDY` should be
kept; they avoid the one-hour timeout failures observed in earlier validation
runs.

## Troubleshooting

### Requested proteins or covariates are missing

Check:

```text
<outdir>/_shared/001_validation/validation_report.txt
```

For protein lists, `--proteins` must be either a comma-separated list of column
names or a file containing one phenotype column name per line.

### Run all proteins with the wrapper

The wrapper default is a 10-protein test list. Use an explicit empty value:

```bash
sbatch src/005_run.sh --proteins ""
```

### Downstream tasks hit time limits

Keep or increase the `time` values for `GWAS_QC`, `META_GROUP`, and `META_STUDY`
in both the base `process` block and `profiles.slurm.process` overrides in
`pipeline/nextflow.config`, then resubmit with `sbatch src/005_run.sh`. The
wrapper uses `-resume`.
