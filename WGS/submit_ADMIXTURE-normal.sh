#!/bin/bash
#SBATCH --job-name=ADMIXTURE
#SBATCH --output=ADMIXTURE-%A.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

source ~/.bashrc
ml load cluster/plink/1.90

## 2.14 Projection Analysis
REF="hgdp_1kg.99AFR_99AMR_99EUR.auto.hg38.maf010.shared.pruned"
STUDY="redlat_tangl.auto.qceed.include_anc.pruned"

## directories
MAINDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25"
DIR="${MAINDIR}/global_ancestry/ADMIXTURE"

cd ${DIR}

start=$(date +%s)

# Verify the two datasets have the same set of SNPs
diff -s "${DIR}/${REF}.bim" "${DIR}/${STUDY}.bim"

# Run unsupervised ADMIXTURE with K=3
admixture "${DIR}/${REF}.bed" 3

# Project entire study sample together
cp "${DIR}/${REF}.3.P" "${DIR}/${STUDY}.3.P.in" # Use learned allele frequencies as (fixed) input
admixture -P "${DIR}/${STUDY}.bed" 3 # Run projection ADMIXTURE 

end=$(date +%s)

## print time elapsed in hours, minutes, and seconds
seconds_elapsed=$((end - start))
hours=$((seconds_elapsed / 3600))
minutes=$(( (seconds_elapsed % 3600) / 60 ))
seconds=$((seconds_elapsed % 60))
echo "Time elapsed: ${hours}h:${minutes}m:${seconds}s"