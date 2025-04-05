# APOE_local_ancestry
Genotyping and local ancestry inference at the APOE locus \
Based off of https://github.com/broadinstitute/gnomad_local_ancestry \
**Input**: phased study genomes (if using array data you will need both the unimputed and imputed/phased versions of your data). For this example, we used array data phased and imputed with the TOPMed server (see https://github.com/HudsonAlpha/Pre-Imputation-QC-Pipeline) \
**Output**: APOE genotypes and local ancestry assignments

## Overview
1. Create a reference panel of AFR, AMR and EUR genomes (PC-AiR)
2. Perform PCA on merged reference panel and study genomes (PC-AiR)
3. Estimate global ancestry (ADMIXTURE)
4. Infer and extract local ancestry of study genomes (RFmix & Tractor)

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

## 4. Infer and extract local ancestry of study genomes (RFmix & Tractor)
- Perform LAI on the imputed/phased study genomes
- Get the haplotype and dosage counts for each ancestry
