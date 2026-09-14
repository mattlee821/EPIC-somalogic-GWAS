# Reproducible genetic loci associated with smoking behaviour

## Scope and interpretation

This report prioritises primary genome-wide association meta-analyses from the Tobacco and Genetics Consortium (TAG), the deCODE/ENGAGE programme, and the GWAS & Sequencing Consortium of Alcohol and Nicotine use (GSCAN), with later large nicotine-dependence meta-analyses used to evaluate replication. “Reproducible” means that the same locus is supported by genome-wide or strong nominal association in an independent primary study/meta-analysis with substantially non-identical samples; the lead SNP need not be identical when different variants are in linkage disequilibrium. This is important because several early studies reported different sentinel SNPs within the same locus.

For direction, I report the effect allele and the direction of its association with the phenotype as analysed in the cited discovery study. For GSCAN binary phenotypes, note that SmkCes was coded so that higher values represented current rather than former smoking.

## Core reproducible loci

| Locus | Gene(s) / interpretation | Lead SNP | Phenotype | Effect direction | Discovery meta-analysis | P-value | Discovery N | Replication / corroboration | Reproducibility |
|---|---|---|---|---|---|---:|---:|---|---|
| 15q25.1 | **CHRNA5–CHRNA3–CHRNB4**; the canonical nicotinic receptor cluster | rs1051730-A | Cigarettes/day | **Increase** (~1.03 CPD/allele) | TAG 2010 | 2.8e-73 | 73,853 | ENGAGE: rs1051730-A, +0.80 CPD, P=2.4e-69, N=76,972; GSCAN: rs16969968-A, P=1.2e-278, N=330,721; iNDiGO: rs16969968, P=1.6e-39 for ND | **Robustly replicated** |
| 15q25.1 | **CHRNA5–CHRNA3–CHRNB4** | rs16969968-A | Nicotine dependence | **Increase in dependence** | iNDiGO 2020, EUR meta-analysis | 1.6e-39 | 46,213 | Prior GWAS and GSCAN/other smoking phenotypes support the same locus; rs16969968 is a functional CHRNA5 D398N variant | **Robustly replicated** |
| 8p11.2 | **CHRNB3–CHRNA6** | rs6474412-T | Cigarettes/day | **Increase** (+0.29 CPD) | ENGAGE/deCODE-led 2010 | 1.4e-08 | 84,956 | Same 8p11 locus was supported for ND in follow-up studies; correlated variants rs1451240/rs13273442 also associated with dependence | **Replicated at locus level** |
| 19q13.2 | **CYP2A6 / CYP2B6** nicotine-metabolism region | rs4105144-C | Cigarettes/day | **Increase** in CPD with C allele in the discovery coding scheme | ENGAGE/deCODE-led 2010 | 2.2e-12 | 83,317 | GSCAN independently replicated the CYP2A6 region; multiple variants, including rs56113850, associated with CigDay (P=4.0e-99). The region also associates with smoking cessation | **Robustly replicated; multiple signals** |
| 19q13.2 | **EGLN2**, adjacent to CYP2A6 | rs3733829-G | Cigarettes/day | **Increase** (+0.33 CPD) | TAG 2010 | 1.0e-08 | 73,853 | Independent follow-up work showed association with CPD and nicotine dependence, and demonstrated that the signal is not simply a proxy for CYP2A6 activity | **Replicated at locus level** |
| 9q34.3 | **DBH** | rs3025343-G | Smoking cessation (former vs current) | **Increase odds of being a former smoker** (OR 1.12) | TAG 2010 | 3.6e-08 | 64,924 | FINRISK, N=11,926, OR 1.12, P=0.094, reproducing the same direction and effect size; later GWAS continue to implicate DBH | **Direction/effect replicated; later GWAS support** |
| 9q34.3 | **DBH** | rs13284520 | Nicotine dependence | **Increase in dependence** for the reported risk allele | iNDiGO 2020, EUR meta-analysis | 1.1e-08 | 46,213 | Same locus also showed very strong CigDay/SmkCes associations in GSCAN | **Replicated at locus level** |
| 20q13.3 | **CHRNA4** | rs151176846 | Nicotine dependence | **Increase in dependence** | iNDiGO 2020, EUR meta-analysis | 1.2e-12 | 46,213 | Independent CHRNA4 GWAS meta-analysis identified rs2273500, P=8.0e-09 after replication; GSCAN smoking traits also implicated CHRNA4 | **Robustly replicated** |
| 20q11.2 | **DNMT3B** | rs910083-C | Nicotine dependence | **Increase in dependence** (OR 1.06 for severe vs mild) | Bierut/Hancock et al. 2017 | 3.7e-08 | 38,602 | Follow-up meta-analysis reproduced the DNMT3B association; nearby **NOL4L** is also implicated in heavy smoking | **Replicated; distinct from NOL4L signal** |
| 20q11.2 | **NOL4L** (near DNMT3B) | rs57342388 indel | Heavy vs never smoking / dependence | **Increase in heaviness/dependence** | UK Biobank heavy-smoking GWAS / subsequent ND meta-analysis | 4.7e-09 in prior heavy-smoking GWAS; P=0.0017 in ND meta-analysis | large UK Biobank study; subsequent ND meta-analysis N=38,602 | Signal remained associated after conditioning on DNMT3B rs910083, indicating two nearby signals | **Replicated at locus level; independent from DNMT3B** |
| 20q13.3 | **CHRNA4** | rs2273500-C | Nicotine dependence | **Increase in dependence** | CHRNA4 GWAS meta-analysis 2016 | 2.3e-06 in discovery; 9.2e-04 in replication; 8.0e-09 combined | 17,074 discovery; 7,469 replication | Same locus also genome-wide significant in iNDiGO; multiple correlated CHRNA4 SNPs | **Robustly replicated** |

### What is especially robust

The most reproducible smoking loci across the prioritised literature are the **15q25 CHRNA5–CHRNA3–CHRNB4 cluster**, **19q13 CYP2A6 region**, **8p11 CHRNB3–CHRNA6**, and **9q34 DBH**. The 20q13 **CHRNA4** and 20q11 **DNMT3B/NOL4L** regions are also reproducible for nicotine dependence/heavy smoking, although their evidence is more phenotype-specific than the 15q25 and 19q13 regions.

## Phenotype-specific findings and important non-replications

### Smoking initiation

TAG 2010 reported eight genome-wide significant SNPs for smoking initiation in the combined analysis; the strongest was the missense **BDNF rs6265-C**, OR 1.06, P=1.8e-08, N=143,023. The same TAG work identified **DBH rs3025343-G** for smoking cessation (OR 1.12, P=3.6e-08, N=64,924). The BDNF result is not as reproducibly established as the canonical 15q25/19q13 loci: later large GSCAN analyses identified hundreds of initiation loci, but individual initiation signals are often cohort- and phenotype-dependent. Consequently, BDNF rs6265 should be treated as a **reported/limited-replication initiation locus**, not as a core cross-consortium smoking locus.

GSCAN’s 2019 meta-analysis had **N=1,232,091** for ever-vs-never regular smoking and identified **378 associated variants**. This very large result established that initiation is highly polygenic. The number of loci is much larger than the stable set seen in pre-GSCAN studies, and many of the new loci have not been independently replicated in another primary meta-analysis of comparable scale.

### Age at smoking initiation

GSCAN analysed **N=341,427** and reported **10 associated variants/loci**. The phenotype is genetically distinct from smoking quantity and cessation, and the paper explicitly reported that the loci were broadly pleiotropic with other smoking traits. However, within the priority TAG/deCODE/iNDiGO literature, I did **not** find a specific age-at-initiation lead SNP with replication strong enough to classify it alongside the canonical reproducible loci above. Therefore, age-at-initiation should presently be regarded as an important **GWAS-discovered but less independently validated** phenotype than CPD, cessation, or composite nicotine dependence.

### Cigarettes per day

The earliest reproducible architecture is especially clear. TAG 2010 identified 15q25, 10q25, and EGLN2-associated signals; ENGAGE/deCODE independently supported 15q25, 8p11, and 19q13. The 19q13 locus is biologically heterogeneous and contains both **CYP2A6/CYP2B6** nicotine-metabolism signals and the nearby **EGLN2** signal. Later GSCAN confirmed the same regions and added many additional loci, including very strong CYP2A6 signals.

### Smoking cessation

The most reproducible early signal is **DBH rs3025343**. TAG found OR=1.12 for former vs current smoking, and an independent FINRISK analysis reproduced the same effect direction and essentially the same effect estimate, although its smaller sample size gave P=0.094. GSCAN later found multiple signals in **CYP2A6**, **DBH**, **CHRNA5/15q25**, and other regions. Thus cessation shares part of its architecture with heaviness/ND but is not genetically identical to them.

### Nicotine dependence

Composite nicotine dependence, typically measured with the Fagerström Test for Nicotine Dependence (FTND), has a smaller but more stable set of loci than smoking initiation. A later synthesis of the primary GWAS literature identified six reproducible genome-wide significant loci: **CHRNB3–CHRNA6 (8p11), DBH (9q34), CHRNA5–CHRNA3–CHRNB4 (15q25), DNMT3B/NOL4L (20q11), and CHRNA4 (20q13)**. The iNDiGO meta-analysis demonstrated especially strong effects at 15q25 and CHRNA4 and corroborated the 9q34 signal.

## Loci reported primarily once or with weak replication

The 2010 TAG/ENGAGE work also nominated several loci that did not become part of the stable replicated set. Examples include **10q25 rs1329650 and rs1028936** for CPD, and **7p14 rs215614/rs215605**. These reached genome-wide significance or near-significance in particular meta-analytic combinations, but their evidence for durable independent replication is substantially weaker than for 15q25, 19q13, 8p11, and DBH. They should therefore be described as **single-study/limited-replication loci**, rather than established smoking loci.

Similarly, many of the hundreds of GSCAN loci—especially those for ever-smoking initiation—are genuine genome-wide discoveries but should not automatically be called “reproducible” merely because they reached P<5e-8 in the GSCAN meta-analysis. GSCAN itself identified **298 genome-wide significant loci** across the four smoking phenotypes, with **259 involving smoking initiation**, illustrating how much of the architecture was discovered de novo in that consortium.

## Summary by phenotype

| Phenotype | Strongest reproducible loci | Overall confidence |
|---|---|---|
| Smoking initiation | Broadly polygenic; limited stable cross-consortium core. BDNF rs6265 is an early reported signal but not as robust as 15q25/19q13. | **Moderate for individual loci; high polygenicity** |
| Age at initiation | GSCAN identified 10 loci; independent replication of individual lead SNPs is limited. | **Lower locus-level replication** |
| Cigarettes/day | **15q25 CHRNA5-A3-B4; 8p11 CHRNB3-A6; 19q13 CYP2A6/CYP2B6; 19q13 EGLN2**; additional GSCAN loci are numerous. | **Very high for the core loci** |
| Smoking cessation | **DBH rs3025343**; also strong overlap with CYP2A6 and 15q25 signals in later GWAS. | **High for DBH; moderate for additional loci** |
| Nicotine dependence | **15q25, 8p11, 9q34, 20q11, 20q13**. | **High for the established core set** |

## Primary sources

1. Tobacco and Genetics Consortium. *Genome-wide meta-analyses identify multiple loci associated with smoking behavior.* Nature Genetics. 2010;42:441–447. DOI: https://doi.org/10.1038/ng.571
2. Thorgeirsson TE et al. *Sequence variants at CHRNB3–CHRNA6 and CYP2A6 affect smoking behavior.* Nature Genetics. 2010;42:448–453. DOI: https://doi.org/10.1038/ng.573
3. Liu M et al. *Association studies of up to 1.2 million individuals yield new insights into the genetic etiology of tobacco and alcohol use.* Nature Genetics. 2019;51:237–244. DOI: https://doi.org/10.1038/s41588-018-0307-5
4. Hancock DB et al. *Genome-wide association study across European and African American ancestries identifies a SNP in DNMT3B contributing to nicotine dependence.* Molecular Psychiatry. 2018. DOI: https://doi.org/10.1038/mp2017-193
5. Hancock DB et al. *Genome-wide meta-analysis reveals common splice site acceptor variant in CHRNA4 associated with nicotine dependence.* [Primary meta-analysis, 2016].
6. iNDiGO Consortium. *Expanding the genetic architecture of nicotine dependence and its shared genetics with multiple traits.* Nature Communications. 2020. DOI: https://doi.org/10.1038/s41467-020-19265-z
7. Bloom AJ et al. *Variants in two adjacent genes, EGLN2 and CYP2A6, influence smoking behavior related to disease risk via different mechanisms.* Human Molecular Genetics. 2014;23:555–561. DOI: https://doi.org/10.1093/hmg/ddt432
8. National FINRISK replication of DBH rs3025343: *Association of the DBH Polymorphism rs3025343 With Smoking Cessation in a Large Population-Based Sample.* Nicotine & Tobacco Research. 2017.

## Key source facts checked

- TAG 2010: n=74,053 overall consortium; CPD meta-analysis showed rs1051730 P=2.8e-73, rs1329650 P=5.7e-10, rs1028936 P=1.3e-9, and EGLN2 rs3733829 P=1.0e-8. Smoking initiation highlighted BDNF rs6265 (P=1.8e-8), while cessation highlighted DBH rs3025343 (P=3.6e-8).
- ENGAGE/deCODE 2010: CPD meta-analysis N=31,266 in the discovery stage, with combined follow-up giving rs6474412 N=84,956, rs4105144 N=83,317, and rs1051730 N=76,972. The 8p11 signal maps to CHRNB3–CHRNA6 and the 19q13 signal to the CYP2A6/CYP2B6 region.
- GSCAN 2019: AgeSmk N=341,427; SmkInit N=1,232,091; CigDay N=337,334; SmkCes N=547,219. Initial GWAS yielded 566 independent genome-wide significant variants in 406 loci across tobacco and alcohol traits; for the four smoking phenotypes, the smoking-specific counts were 10, 378, 55, and 24 associated variants respectively.
- iNDiGO: European-ancestry ND meta-analysis N=46,213; strongest established signals included rs16969968 at CHRNA5 (P=1.6e-39), rs151176846 at CHRNA4 (P=1.2e-12), and rs13284520 at DBH (P=1.1e-8).

## Bottom line

For applications such as Mendelian randomization, genetic-risk scoring, or adjustment for smoking behaviour, the safest “established” smoking loci are **CHRNA5/CHRNA3/CHRNB4 (15q25), CYP2A6/CYP2B6 and EGLN2 (19q13), CHRNB3/CHRNA6 (8p11), and DBH (9q34)**, with **CHRNA4 (20q13)** and **DNMT3B/NOL4L (20q11)** added when the target phenotype is nicotine dependence/heavy smoking. Smoking initiation and especially age at initiation are substantially more polygenic, and most GSCAN-discovered loci should be regarded as discovery loci unless independently replicated.
