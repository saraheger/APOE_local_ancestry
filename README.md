# APOE_local_ancestry
Genotyping and local ancestry inference at the APOE locus \
Based off of https://github.com/broadinstitute/gnomad_local_ancestry
## Overview
1. Create a reference panel of AFR, AMR and EUR genomes (PC-AiR)
2. Perform PCA on merged reference panel and study genomes (PC-AiR)
4. Estimate global ancestry (ADMIXTURE)
5. Impute & phase study genomes (TOPMed & Eagle)
6. Infer local ancestry of study genomes (RFmix)
7. Extract allele ancestries (Tractor)

## 1. Create a reference panel of AFR, AMR and EUR genomes (PC-AiR)
- Download phased reference genomes BCF files: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
- Convert to plink & keep only AFR, AMR and EUR population genomes
- Run PC-AiR and select similar numbers of genomes per ancestral population
