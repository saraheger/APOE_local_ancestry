#!/bin/bash
#SBATCH --job-name=tractor
#SBATCH --output=tractor-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

ml load bcftools
ml load python

# get chromosome
CHR=${SLURM_ARRAY_TASK_ID}
echo "chr${CHR}"

## directories
MAIN_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25"
INDIR1="${MAIN_DIR}/SHAPEIT5_v2/3_phase_common"
INDIR2="${MAIN_DIR}/LAI/RFmix_v2"
OUTDIR1="${MAIN_DIR}/LAI/RFmix_v2/vcfs"
OUTDIR2="${MAIN_DIR}/LAI/RFmix_v2/Tractor"
mkdir -p ${OUTDIR1}
mkdir -p ${OUTDIR2}
cd ${OUTDIR2}

## files
QUERY_BCF="${INDIR1}/redlat_tangl.chr${CHR}.qceed.include_anc.maf0001.shapeit5_common.bcf"
QUERY_VCF="${OUTDIR1}/redlat_tangl.chr${CHR}.qceed.include_anc.maf0001.shapeit5_common.vcf.gz"
MSP="${INDIR2}/redlat_tangl.chr${CHR}.qceed.include_anc.maf0001.shapeit5_common.LAI.msp.tsv"

start=$(date +%s)

# BCF to VCF
bcftools view ${QUERY_BCF} \
    --write-index \
    -Oz -o ${QUERY_VCF}

# extract tracts
python /cluster/home/seger/software/Tractor/scripts/extract_tracts1.py \
    --msp ${MSP} \
    --vcf ${QUERY_VCF} \
    --num-ancs 3 \
    --compress-output \
    --output-dir ${OUTDIR2}

end=$(date +%s)

## print time elapsed in hours, minutes, and seconds
seconds_elapsed=$((end - start))
hours=$((seconds_elapsed / 3600))
minutes=$(( (seconds_elapsed % 3600) / 60 ))
seconds=$((seconds_elapsed % 60))
echo "Time elapsed: ${hours}h:${minutes}m:${seconds}s"
# chr19 = 5h:3m:26s