# APOE_local_ancestry
Genotyping and local ancestry inference at the APOE locus
Based off of https://github.com/broadinstitute/gnomad_local_ancestry
## Overview
1. Create a reference panel of AFR, AMR and EUR genomes using PC-AiR
2. Run PC-AiR on reference panel genomes and study genomes
3. Run ADMIXTURE projection analysis for global ancestry
4. Phase (and impute if array data) study genomes
5. Run RFmix on phased study genomes
6. Extract allele ancestries
