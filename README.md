# APOE_local_ancestry
Genotyping and local ancestry inference at the APOE locus \
Based off of https://github.com/broadinstitute/gnomad_local_ancestry
## Overview
1. Create a reference panel of AFR, AMR and EUR genomes (PC-AiR)
2. Perform PCA on merged reference panel and study genomes (PC-AiR)
3. Estimate global ancestry (ADMIXTURE)
4. Impute & phase study genomes (TOPMed & Eagle)
5. Infer local ancestry of study genomes (RFmix)
6. Extract allele ancestries (Tractor)

## 1. Create a reference panel of AFR, AMR and EUR genomes (PC-AiR)
- Download phased reference genomes BCF files: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
- Convert to plink & keep only AFR, AMR and EUR population genomes
- Run PC-AiR and select similar numbers of genomes per ancestral population

## 2. Perform PCA on merged reference panel and study genomes (PC-AiR)
- Find the variants shared between the reference genomes and the unimputed study genomes
- merge with PLINK
- Run PC-AiR to visualize where the admixed population clusters 

## 3. Estimate global ancestry (ADMIXTURE)
- Prune the reference panel genomes (including only the variants shared with the study genomes)
- Extract the same prune set from the study genomes
- Run ADMIXTURE on the reference panel genomes
- Project the study genomes onto the reference panel output

## 4. Impute & phase study genomes (TOPMed & Eagle)
- Used this pipeline for array data: https://github.com/HudsonAlpha/Pre-Imputation-QC-Pipeline

## 5. Infer local ancestry of study genomes (RFmix)
- Perform LAI on the imputed/phased study genomes

## 6. Extract allele ancestries (Tractor)
- Get the haplotype and dosage counts for each ancestry
