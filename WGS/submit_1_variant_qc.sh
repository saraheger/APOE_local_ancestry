#!/bin/bash
#SBATCH --job-name=1_variant_qc
#SBATCH --output=1_variant_qc-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=30
#SBATCH --array=1-2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

ml load cluster/bcftools/1.18

# get threads
threads=${SLURM_CPUS_PER_TASK}

# get chromosome
CHR=${SLURM_ARRAY_TASK_ID}
echo "chr${CHR}"

## directories
INDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/chrs_vcfs"
OUTDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5/1_variant_qc"
mkdir -p "${OUTDIR}"

## files
IN="${INDIR}/redlat_tangl.chr${CHR}.vcf.gz"
OUT="${OUTDIR}/redlat_tangl.chr${CHR}.qceed.bcf"
OUT2="${OUTDIR}/redlat_tangl.chr${CHR}.NOqceed.stats.bcf"
TMP="${OUTDIR}/redlat_tangl.chr${CHR}.tmp.bcf"

# Splits multiallelic variants into multiple biallelic records
# Removes the non-GT items
# Computes Fraction of missing genotypes
# Keep PASS & remove symbolic alleles, missing > 5%, and monomorphic variants

bcftools norm --threads ${threads} -m -any ${IN} -Ou | \
bcftools annotate -x "FORMAT" -Ou | \
bcftools +fill-tags --threads ${threads} -Ob -o ${TMP} -- -t F_MISSING && \
bcftools view ${TMP} --threads ${threads} -f PASS -e 'ALT="*" | F_MISSING>0.05 | INFO/AF=0 | INFO/AF=1' -Ob -o ${OUT} && \
bcftools index ${OUT} --threads ${threads} && \
bcftools view -G ${TMP} --threads ${threads} -Ob -o ${OUT2} && \
bcftools index ${OUT2} --threads ${threads} && \
rm ${TMP}


