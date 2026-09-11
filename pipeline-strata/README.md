# EPIC SomaScan protein GWAS — stratified extended analysis

Protein GWAS run **within strata of a participant-level variable** (ever vs never
smokers, BMI < 30 vs ≥ 30, and anything else you define), adjusting for disease
status as a covariate.

This is a **fork, not a branch**, of the publication pipeline. The published
analysis is frozen at git tag **`paper-v1`** and nothing in this tree modifies
it. If you need to reproduce the paper, check out that tag and use `pipeline/`
and `src/` at the repository root — not this directory.

**Smoking is verified and ready.** `smoke_statVV` → `SMOKE_CAT`:

| raw | n | → | stratum |
|---|---|---|---|
| `Never` | 8,046 | → | `NEVER` |
| `Former` | 4,911 | → | `EVER` |
| `Smoker` | 4,535 | → | `EVER` |
| `Unknown` | 102 | → | *excluded from both* |

(Full-cohort counts; the genotyped analysis set is smaller after the
proteomics intersect and cross-study dedup.)

**BMI is declared but off.** Its source column name is not yet confirmed, so
`BMI_CAT` is marked `optional` — 003 skips it with a warning and the
`bmi_lo`/`bmi_hi` strata simply never run. Add the real column name to
`strata_spec$BMI_CAT$source` in `src/003_data-phenotype.R` to switch it on;
nothing else needs changing.

---

## What a "stratum" is here

A stratum is a **named whitelist of accepted values in one column**, declared in
`pipeline/params.yaml`:

```yaml
strata:
  all:    { column: "",          values: [] }              # everyone
  ever:   { column: "SMOKE_CAT", values: ["EVER"] }
  never:  { column: "SMOKE_CAT", values: ["NEVER"] }
  bmi_lo: { column: "BMI_CAT",   values: ["LO"] }
  bmi_hi: { column: "BMI_CAT",   values: ["HI"] }
pca_group: "all"
```

Membership is **always explicit**. A value that no stratum lists — `Unknown`,
a blank, a level you forgot — belongs to **no** stratum rather than falling into
a complement group.

This is the substantive difference from the frozen pipeline, which defined its
control arm as `!= cases_value`. That is safe for a binary disease code and
wrong for anything else. On the real smoking data, `!= Never` gives an ever arm
of **9,548** against a true **9,446** — the 102 `Unknown` silently analysed as
smokers.
The frozen code also cannot run this analysis at all: it resolves the grouping
column *after* covariate subsetting has already dropped it, so a strata column
living in the phenofile raises `Group column not found`.

## Adjusting for disease status

`PHENO` exists **only in the study `.psam`** — the frozen `003` deliberately
keeps it out of the phenofile so it cannot become a trait. The frozen
`prepare_covariates.py` pulls only the *grouping* column from the `.psam`, so
`PHENO` was reachable as a covariate only while it was itself the grouping
variable. Here the grouping variable is smoking or BMI, so `prepare_covariates.py`
resolves **any** requested covariate against the `.psam` as a fallback. Request
it explicitly:

```
--covariates "Sex,Age_Recr,PHENO"
```

which is the default in `src/005_run.sh`. A requested covariate that exists in
neither file is a hard error, not a silent drop.

## Ancestry PCs

PCs are computed once per study on `pca_group` — the **unfiltered** stratum —
and shared with every other stratum. They are never computed inside a stratum:
doing so within an exposure-defined subgroup can absorb real association signal
and makes strata non-comparable. The pipeline refuses to start if `pca_group`
names a filtered stratum, and forces `pca_group` into every study's group list
(omitting it would produce no PCs, join empty at `ADD_PCS`, and silently emit
zero results).

## What is shared with the frozen tree, and what is not

**Read, shared** — re-deriving these would cost days for byte-identical output,
so this tree points at them and never writes to them:

| Artefact | Produced by |
|---|---|
| Per-study PLINK2 filesets (`GWAS_GENETICS_DIR`) | root `src/002_data-genetics.sh` |
| Conda environments (`env/`) | root `src/001_source.sh` |
| LDSC reference data, protein–gene map | root `src/004_ldsc-herit.sh` |
| SomaScan RDS (`GWAS_EXTERNAL_SOMA_RDS`) | upstream |

**Write, separate** — every path this tree writes to is distinct, and several
guards enforce it at runtime:

| Variable | Note |
|---|---|
| `GWAS_PHENO_DIR` | **New.** Explicit, never derived. The frozen `003` derives its phenotype dir from `GWAS_GENETICS_DIR`; since both trees share that directory, deriving it here would overwrite the published phenofile. `003` aborts if it resolves to the frozen path. |
| `GWAS_RESULTS_DIR`, `GWAS_WORK_DIR`, `GWAS_STAGING_DIR`, `GWAS_SAMPLESHEET` | `-strata` suffixed. `005_run.sh` aborts if `GWAS_WORK_DIR` is the frozen tree's. |

Two variables are **new and required**, because the frozen modules resolved
these relative to their own tree (`${projectDir}/../env`, `${projectDir}/tools`)
— which here would point at `pipeline-strata/env` and `pipeline-strata/pipeline/tools`,
neither of which exists, and every process would fail to start:

| Variable | Points at |
|---|---|
| `GWAS_ENV_DIR` | the frozen tree's `env/` (conda environments from `001`) |
| `GWAS_METAL_BIN` | the frozen tree's `pipeline/tools/metal/metal` |



## Running — this analysis starts at 003

Everything upstream of 003 already exists and is reused: the per-study PLINK2
filesets, the conda environments, the METAL binary, and the LDSC reference data.
`001_source.sh`, `002_data-genetics.sh` and `004_ldsc-herit.sh` are present for
completeness but **refuse to run** from this tree (override with
`GWAS_ALLOW_UPSTREAM_RERUN=1` only if you are certain). Re-running `002` in
particular would overwrite genetics filesets that the frozen analysis reads.

```bash
cp .env.example .env          # then check every path
sbatch  src/003_data-phenotype.sh    # derives SMOKE_CAT + writes the samplesheet
sbatch  src/005_run.sh --out /WHERE/epic-gwas/strata-test \
                       --proteins "" --studies "Neuro_01,Brea_01_Erneg"
bash    src/006_staging.sh --out <outdir> --transfer <dir>
```

Read 003's log before launching 005: it prints the raw-to-mapped crosstab for
every derived variable and the per-study stratum counts that decide which
strata actually run.

`003` decides which strata each study is large enough to run (`min_stratum_n`,
read back out of `params.yaml` so the samplesheet can never disagree with the
pipeline) and writes them pipe-separated into `analysis_groups`.

Useful `005_run.sh` flags: `--covariates`, `--meta-study` (strata to
meta-analyse across studies, each independently), `--meta-group` (within-study
meta across a *pair* of strata — rarely wanted, since pooling ever+never just
reconstructs `all`), `--samplesheet`.

## Output layout

```
<outdir>/meta/<feature>/<stratum>/013_QC/meta.tsv.gz
<outdir>/meta/<feature>/<stratum>/014_metrics/        # LDSC + heritability
<outdir>/all_metrics.tsv                              # rows carry study + group
```

Staged for transfer as `<transfer>/<stratum>/<feature>.tsv.gz`.

The stratum had to be added to these publish paths: the frozen modules wrote to
`meta/<feature>/013_QC/` with `overwrite: true`, which was safe only because a
single group was ever meta-analysed. Meta-analysing several strata in one run
would have had them silently clobber each other, last-to-finish winning.

## Two caveats on interpretation

**Disease adjustment may be partly redundant.** The input proteins come from
`NormalizedSoma_PlateCorr*AccountforDisease*_InvRank` — disease status is
already regressed out during normalisation. Adjusting for `PHENO` on top is
defensible but worth checking against what that normalisation actually did.

**Inverse-rank normalisation was applied to the full sample, not within
strata.** Within ever/never the phenotype is therefore not exactly standard
normal. Harmless for REGENIE's linear model, and it keeps effect sizes
comparable across strata — but do not describe them as within-stratum SD units.
