# Genomic Loci and Implicated Genes Reproducibly Associated with Smoking Behaviors

This report synthesizes findings from major large-scale genome-wide association studies (GWAS) and meta-analyses conducted by international consortia, including the **GWAS & Sequencing Consortium of Alcohol and Nicotine use (GSCAN)**, the **Tobacco and Genetics (TAG) Consortium**, and **deCODE genetics**. 

It outlines primary genomic loci associated with five core smoking phenotypes:
1. **Smoking Initiation (Init):** Ever vs. never regular smoker.
2. **Age of Smoking Initiation (AgeInit):** Age at which regular smoking commenced.
3. **Cigarettes Per Day (CPD):** Average daily cigarette consumption among current/former smokers.
4. **Smoking Cessation (Cess):** Former vs. current smokers (successful cessation).
5. **Nicotine Dependence (ND):** Measured via Fagerström Test for Nicotine Dependence (FTND) or DSM criteria.

---

## 1. Replicated Loci Across Independent Large-Scale Meta-Analyses

The following loci have met genome-wide significance ($P < 5 \times 10^{-8}$) and have been independently validated across multiple primary consortia (deCODE, TAG, GSCAN Phase 1 & 2, UK Biobank, and Million Veteran Program):

| Cytoband | Implicated Gene(s) | Lead SNP (rsID) | Trait / Phenotype | Risk/Effect Allele & Direction | Discovery $P$-Value | Sample Size ($N$) | Primary Study / Consortium | Replication Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **15q25.1** | *CHRNA5-CHRNA3-CHRNB4* | **rs16969968** (or rs1051730) | CPD, ND, Cess | **A allele (D398N):** Increases CPD & ND severity; decreases cessation success | $1.1 \times 10^{-247}$ | ~263,954 | GSCAN Phase 1 (Liu 2019) / deCODE (Thorgeirsson 2008) | **Replicated** (deCODE, TAG, GSCAN 1/2, MVP) |
| **19q13.2** | *CYP2A6 / EPHX2* | **rs56113850** (or rs4105144) | CPD, ND | **C allele:** Associated with faster nicotine metabolism & increased daily CPD | $2.0 \times 10^{-42}$ | ~263,954 | GSCAN Phase 1 (Liu 2019) / TAG (Furberg 2010) | **Replicated** (TAG, GSCAN 1/2, MVP) |
| **8p11.21** | *CHRNB3-CHRNA6* | **rs6474412** | CPD, Init | **A allele:** Increases CPD and odds of regular smoking initiation | $5.3 \times 10^{-23}$ | ~263,954 | GSCAN Phase 1 (Liu 2019) / TAG (Furberg 2010) | **Replicated** (TAG, GSCAN 1/2) |
| **3p21.1** | *CADM2* | **rs62250714** | Init | **A allele:** Increases risk of smoking initiation (risk-taking behavior) | $1.3 \times 10^{-44}$ | ~1,232,091 | GSCAN Phase 1 (Liu 2019) / UK Biobank | **Replicated** (UK Biobank, GSCAN 1/2, MVP) |
| **11q23.2** | *NCAM1-TTC12-ANKKK1-DRD2* | **rs11214607** | Init, CPD | **A allele:** Increases initiation likelihood and heavier smoking | $2.4 \times 10^{-20}$ | ~1,232,091 | GSCAN Phase 1 (Liu 2019) / TAG (2010) | **Replicated** (TAG, GSCAN 1/2) |
| **9q13.2** | *DBH* | **rs3025343** | Cess, CPD | **A allele:** Reduces probability of cessation (higher persistence) | $3.6 \times 10^{-12}$ | ~547,219 | GSCAN Phase 1 (Liu 2019) / TAG (2010) | **Replicated** (TAG, GSCAN 1/2) |
| **7p14.3** | *AHR* | **rs4410790** | AgeInit, CPD | **C allele:** Associated with younger age at smoking initiation | $2.1 \times 10^{-11}$ | ~341,427 | GSCAN Phase 1 (Liu 2019) / TAG (2010) | **Replicated** (TAG, GSCAN 1/2) |
| **11p14.1** | *BDNF* | **rs6265** (Val66Met) | Init, AgeInit | **T allele:** Increases initiation propensity & affects onset age | $4.1 \times 10^{-14}$ | ~1,232,091 | GSCAN Phase 1 (Liu 2019) / TAG (2010) | **Replicated** (TAG, GSCAN 1/2) |

---

## 2. Loci Reported Solely in Specific Cohorts or Single Discovery Analyses

While major nicotinic receptor and pharmacokinetic loci are robustly replicated across multi-ancestry cohorts, several loci are context-specific, reported in single discovery meta-analyses, or exhibit population-specific architecture:

* **11q13.2 (*SLC25A21*) — rs77864389**
  * **Phenotype:** Nicotine Dependence (ND) / CPD
  * **Lead SNP:** rs77864389 ($P = 4.2 \times 10^{-9}, N \approx 29,000$)
  * **Consortium/Study:** GSCAN Phase 2 African-ancestry discovery cohort (Saunders et al., 2022).
  * **Replication Status:** **Reported once / Population-specific.** Allele frequency is near zero in European populations ($MAF < 0.001$), failing to reach significance in European-majority meta-analyses.

* **1q21.3 (*PDE4B*) — rs7598889**
  * **Phenotype:** Nicotine Dependence severity
  * **Lead SNP:** rs7598889 ($P = 2.8 \times 10^{-8}, N \approx 148,000$)
  * **Consortium/Study:** COGEND / Early TAG subgroup candidate GWAS.
  * **Replication Status:** **Reported once / Unconfirmed.** Did not clear conditional genome-wide significance thresholds in larger GSCAN Phase 1 or 2 multi-ancestry meta-analyses.

* **20q13.33 (*GSS*) — rs6071980**
  * **Phenotype:** Smoking Cessation (Cess)
  * **Lead SNP:** rs6071980 ($P = 3.1 \times 10^{-8}, N = 74,053$)
  * **Consortium/Study:** TAG Consortium (Furberg et al., 2010).
  * **Replication Status:** **Unreplicated.** Marginal association that failed to meet $P < 5 \times 10^{-8}$ in subsequent GSCAN Phase 1 ($N = 547,219$) updates.

---

## 3. Biological Pathways and Functional Mechanisms

The GWAS loci partition into distinct functional domains according to smoking behavioral stage:

1. **Direct Nicotinic Acetylcholine Receptor Signaling (*CHRNA5-CHRNA3-CHRNB4*, *CHRNB3-CHRNA6*):**
   * Encodes $\alpha5$, $\alpha3$, $\beta4$, $\beta3$, and $\alpha6$ subunits of neuronal nicotinic acetylcholine receptors (nAChRs).
   * The missense mutation **rs16969968** (D398N) in *CHRNA5* alters receptor desensitization kinetics in the habenulo-interpeduncular pathway, blunting nicotine aversion and directly driving heavy consumption (CPD) and physiological dependence.

2. **Pharmacokinetics & Nicotine Metabolism (*CYP2A6*):**
   * Encodes Cytochrome P450 2A6, the primary hepatic enzyme responsible for metabolizing nicotine into cotinine.
   * Genetic variants leading to fast nicotine metabolism drive individuals to smoke more cigarettes per day to maintain steady-state circulating nicotine levels.

3. **Dopaminergic & Neurodevelopmental Circuitry (*DRD2/ANKKK1*, *CADM2*, *BDNF*, *DBH*):**
   * *CADM2* (Cell Adhesion Molecule 2) and *BDNF* (Brain-Derived Neurotrophic Factor) mediate executive function, behavioral disinhibition, and impulsivity, strongly predisposing individuals to **Smoking Initiation**.
   * *DRD2/ANKK1* regulates dopamine receptor availability in reward pathways, influencing both early initiation and long-term nicotine reinforcement.
   * *DBH* (Dopamine Beta-Hydroxylase) regulates the conversion of dopamine to norepinephrine, impacting the neurochemical dynamics of **Smoking Cessation**.

---

## Key References & Primary Consortia
- **GSCAN Phase 1:** Liu et al. (2019). *Association analyses of up to 1.2 million individuals yield new insights into the genetic etiology of tobacco and alcohol use.* Nature Genetics, 51(2), 237–244.
- **GSCAN Phase 2:** Saunders et al. (2022). *Human genetics of tobacco and alcohol use in multi-ancestry populations.* Nature, 612(7941), 720–728.
- **TAG Consortium:** Furberg et al. (2010). *Genome-wide association meta-analyses identify multiple loci associated with smoking behavior.* Nature Genetics, 42(5), 441–447.
- **deCODE Genetics:** Thorgeirsson et al. (2008). *A variant associated with nicotine dependence, lung cancer and peripheral arterial disease.* Nature, 452(7187), 638–642.
