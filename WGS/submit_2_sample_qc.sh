#!/bin/bash
#SBATCH --job-name=2_sample_qc
#SBATCH --output=2_sample_qc-%A_%a.out
#SBATCH --time=4:00:00
#SBATCH --partition=normal
#SBATCH --mem=20G
#SBATCH --cpus-per-task=30
#SBATCH --array=22
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

## SHAPEIT5_v2 changes:
# using only ['AMR', 'CAH', 'EUR', 'AAC', 'AFR'] samples (/cluster/home/seger/scripts/redlat_tangl/global_ancestry/genotools_include_anc.py)
# using only trios & duos that were confirmed by KING/PC-Relate (on kosik)

ml load bcftools

# get threads
threads=${SLURM_CPUS_PER_TASK}

# get chromosome
CHR=${SLURM_ARRAY_TASK_ID}
echo "chr${CHR}"

## directories
ADMIXTURE_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/global_ancestry/ADMIXTURE"
INDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5_v3/1_variant_qc"
OUTDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5_v3/2_sample_qc"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

## files
KEEPSMP="${ADMIXTURE_DIR}/include_anc_sampleIDs.txt"

IN="${INDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.bcf"
OUT="${OUTDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.include_anc.bcf"

bcftools view --threads ${threads} -S ${KEEPSMP} ${IN} -Ob -o ${OUT} && \
bcftools index --threads ${threads} ${OUT}
