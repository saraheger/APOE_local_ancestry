#!/bin/bash

# reference panel n = 309

## 2.14 Projection Analysis
REF="hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled.auto.hg38.ext.pruned"
STUDY="E280A_340.auto.hg38.ext.pruned"
RELS="E280A_rels.auto.hg38.ext.pruned"
UNRELS="E280A_unrels.auto.hg38.ext.pruned"

## directories
DIR="/scratch/eger/projects/E280A_SNP_arrays/E280A_340/ADMIXTURE_v1"

cd ${DIR}

start=$(date +%s)

# Verify the two datasets have the same set of SNPs
diff -s "${DIR}/${REF}.bim" "${DIR}/${STUDY}.bim"

# Run unsupervised ADMIXTURE with K=3
admixture "${DIR}/${REF}.bed" 3

# Converged in 12 iterations (65.917 sec)
# Loglikelihood: -22532337.621480
# Fst divergences between estimated populations: 
#         Pop0    Pop1
# Pop0
# Pop1    0.099
# Pop2    0.129   0.091
# Writing output files.

## Version 1:
# Project entire study sample together
cp "${DIR}/${REF}.3.P" "${DIR}/${STUDY}.3.P.in" # Use learned allele frequencies as (fixed) input
admixture -P "${DIR}/${STUDY}.bed" 3 # Run projection ADMIXTURE 

## Version 2:
# Project the Distantly Related Group First
cp "${DIR}/${REF}.3.P" "${DIR}/${UNRELS}.3.P.in" # Use learned allele frequencies as (fixed) input
admixture -P "${DIR}/${UNRELS}.bed" 3 # Run projection ADMIXTURE 

# Project the Very Related Group Using the Updated Model
cp "${DIR}/${UNRELS}.3.P" "${DIR}/${RELS}.3.P.in"
admixture -P "${DIR}/${RELS}.bed" 3

end=$(date +%s)

## print time elapsed in hours, minutes, and seconds
seconds_elapsed=$((end - start))
hours=$((seconds_elapsed / 3600))
minutes=$(( (seconds_elapsed % 3600) / 60 ))
seconds=$((seconds_elapsed % 60))
echo "Time elapsed: ${hours}h:${minutes}m:${seconds}s"
