#!/bin/bash
#SBATCH --job-name=1_variant_qc
#SBATCH --output=1b_variant_qc-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --array=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

ml load bcftools
ml load cluster/plink/1.90

# get threads
threads=${SLURM_CPUS_PER_TASK}

# get chromosome
CHR=${SLURM_ARRAY_TASK_ID}
echo "chr${CHR}"

## directories
PED_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/GENESIS_v2"
INDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5/1_variant_qc"
OUTDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5_v3/1_variant_qc"
mkdir -p "${OUTDIR}"
cd ${OUTDIR}

## files
PED="${PED_DIR}/redlat_tangl_confirmed.bcftools.ped"
IN="${INDIR}/redlat_tangl.chr${CHR}.qceed.bcf"
OUT="${OUTDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.bcf"

# skipping: GQ > 20 for now...
# only SNPs
# set mendelian errors as missing
# bcftools view --threads ${threads} ${IN} -v snps -Ou | \
# bcftools +mendelian2 -P ${PED} -m ad -Ob -o ${OUT} && \
# bcftools index --threads ${threads} ${OUT}

plink --bcf ${OUT} \
    --keep-allele-order \
    --const-fid \
    --allow-no-sex \
    --make-bed \
    --out redlat_tangl.chr${CHR}.qceed.snps.mend