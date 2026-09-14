# 008 — comparing pQTL effects between strata

Answers the question the pipeline sets up but does not itself test: **does a
pQTL act differently in one stratum than another?**

Reads the staged tree from `src/006_staging.sh`:

```
<staged>/<seqid>/{ever,never,meta}.tsv.gz
```

```bash
bash src/008_comparison/run_comparison.sh \
  --staged /WHERE/epic-gwas/transfer-strata \
  --out    analysis/comparison
```

Reference data (LD panel, gene map, HLA window, cis window, clump r²) is read
out of `pipeline/params.yaml`, so it cannot drift from what the pipeline used.

---

## Method

### 1. Discovery comes from the reference, never from a stratum

Selecting pQTLs in `ever` and then comparing them to `never` produces **winner's
curse**: variants chosen for a large effect regress toward the mean in the other
arm, manufacturing differences out of noise. This is the standard way to
generate false "exposure-specific" pQTLs.

Selecting on an inverse-variance-weighted meta is *exactly* unbiased for the
difference. With sampling variances `V` and `W = 1/Vₐ + 1/V_b`:

```
Cov(β_meta, βₐ − β_b) = (1/W)·[ Vₐ/Vₐ − V_b/V_b ] = 0
```

identically, for any variances, because IVW weights are ∝ 1/V. METAL runs with
`SCHEME STDERR` (i.e. IVW), and fixed-effect IVW is associative, so the
two-stage construction — within study across strata, then across studies —
equals one-stage over all study × stratum estimates.

**Blind spot.** A pQTL present in only one stratum has its effect roughly halved
in the meta and may miss the threshold — exactly the case of interest. Those
variants are reported as a **secondary** tier, flagged, because their difference
statistics are selection-biased. Never quote a secondary p-value as a test.

### 2. Independence by LD clumping

Sumstats carry ~11M correlated variants; treating them as independent tests is
meaningless. Clumped with PLINK against the same LD reference the heritability
module uses, HLA dropped (extended LD makes clumping unreliable there). Falls
back to distance-based greedy pruning when no panel is available.

### 3. The comparison

```
δ = βₐ − β_b        SE(δ) = √(SEₐ² + SE_b²)        z = δ / SE(δ)
```

Valid because the strata are **disjoint individuals**, so the estimates are
independent and variances add (Altman & Bland's interaction test).

This is *stronger* than a pooled `SNP × stratum` term, which assumes common
covariate effects and a common residual variance across strata. Here each
stratum is fitted entirely separately, so neither is assumed.

### 4. Global test first, per-locus second

At ~4,000 per stratum, `SE(δ) = √2 × SE`, so the difference test has about half
the power of the main effect:

| MAF | detectable δ (p<0.05) | (Bonferroni, 100 tests) |
|---|---|---|
| 0.1 | 0.103 | 0.183 |
| 0.2 | 0.077 | 0.137 |
| 0.5 | 0.062 | 0.110 |

A strong cis-pQTL is 0.3–0.5 SD, so this detects **near-absence in one arm**,
not subtle modulation. Asking whether the *whole* z distribution departs from
N(0,1) is far better powered than any single locus, so that is reported first.
Random control proteins (`--annot`, `source=random_control`) provide the
empirical null: under H₀, `var(z) = 1`.

### 5. Two QC gates

- **EAF concordance** — the strata are the same populations, so allele
  frequencies must agree. A discrepancy signals a strand or allele problem that
  would otherwise surface as a spurious effect difference.
- **Across-study I² within a stratum** — high heterogeneity means that stratum's
  estimate is unstable, and a difference against the other may be one discordant
  cohort rather than the exposure.

Allele harmonisation is applied before any differencing: `ever` and `never` are
separate METAL runs, so the same variant can return with `EA`/`OA` swapped.
Alleles are upper-cased (METAL lower-cases them), matched as a set, `BETA`
negated and `EAF` complemented where flipped, and anything unmatched dropped.

---

## Outputs

| file | contents |
|---|---|
| `pqtl_comparison.tsv` | one row per index variant: per-stratum β/SE/P/N, δ, z, p, Bonferroni p, cis/trans, QC flags, tier |
| `global_test.tsv` | n, mean z, var z, KS p, λ and the p<0.05 fraction for targeted / control / secondary sets |
| `comparison_plots.png` | β-vs-β with the line of equality; QQ of difference p-values with controls overlaid; z against N(0,1) |
| `comparison.log` | PLINK clumping output |

Multiple-testing correction is applied over **primary, QC-passing,
non-control** variants only — those are the actual hypothesis tests.

## Validation

Checked against synthetic data with known truth:

- A planted stratum-specific effect (β=0.40 vs 0.00) is recovered at
  δ=0.36, Bonferroni p=9e-17; a planted equal effect (0.35 vs 0.34) is not
  flagged.
- Under a pure null of 40 proteins × 6 shared pQTLs: n=235, mean z=−0.013,
  **var z=1.00**, KS p=0.84, exactly 5% below p<0.05, zero false positives.
- Both fixtures allele-flip 30% of `never` variants. Broken harmonisation
  inflates `var(z)` dramatically, so these numbers also validate that step.

## Running it on the whole cluster

`run_comparison.sh` runs a multiprocessing pool, which is capped at one node's
cores. For a full panel use the array submitter instead, from the login node:

```bash
bash src/008_comparison/submit_comparison.sh \
  --annot /data/.../phenotype-strata/smoking-panel.annot.tsv
```

It submits two jobs:

| phase | shape | why |
|---|---|---|
| map | `--array=0-N`, 1 cpu each | proteins are independent; wall time is the slowest protein, not the sum |
| gather | one job, `afterok` on the array | Bonferroni denominators and the var(z)/F tests need every protein in one frame |

Each map task writes `<out>/shards/part_NNNNN.pkl`; the gather job combines them
and writes the final tables. Shards are assigned by STRIDE (`proteins[k::n]`),
not contiguous blocks, because proteins differ several-fold in how many variants
clear the threshold and striding spreads the heavy ones across tasks.

`--gather` refuses to run if any shard is missing. A silently absent shard drops
proteins from the global tests, and nothing downstream could detect it.

Useful flags: `--per-task N` (proteins per task, default 1), `--throttle N` (cap
concurrent tasks), `--dry-run`.
