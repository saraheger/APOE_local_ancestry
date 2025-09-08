#!/bin/bash
#SBATCH --job-name=rfmix
#SBATCH --output=rfmix-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --array=14
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

# chr1: highmem, 800G, 80 CPU
# chr20: normal, 80G, 40 CPU

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
RESOURCE_DIR="/cluster/home/seger/references/gnomad/resources"
MAIN_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25"
INDIR="${MAIN_DIR}/SHAPEIT5_v2/3_phase_common"
OUTDIR="${MAIN_DIR}/LAI/RFmix_v2"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

## files
GENETIC_MAP="${RESOURCE_DIR}/1kg/1000G_genetic_map/genetic_map_hg38_withX_map_only.txt"
SAMPLE_MAP="${RESOURCE_DIR}/hgdp_1kg/ADMIXTURE_99AFR_99AMR_99EUR/hgdp_1kg.99AFR_99AMR_99EUR_samplemap.txt"
REF_BCF="${RESOURCE_DIR}/hgdp_1kg/hgdp1kgp_chr${CHR}.filtered.SNV_INDEL.phased.shapeit5.bcf"
QUERY_BCF="${INDIR}/redlat_tangl.chr${CHR}.qceed.include_anc.maf0001.shapeit5_common.bcf"
OUT="redlat_tangl.chr${CHR}.qceed.include_anc.maf0001.shapeit5_common.LAI"

start=$(date +%s)

rfmix \
    -f ${QUERY_BCF} \
    -r ${REF_BCF} \
    -m ${SAMPLE_MAP} \
    -g ${GENETIC_MAP} \
    -o ${OUT} \
    -e 2 \
    -G 20 \
    --reanalyze-reference \
    --chromosome=chr${CHR} \
    --n-threads=${threads}

end=$(date +%s)

## print time elapsed in hours, minutes, and seconds
seconds_elapsed=$((end - start))
hours=$((seconds_elapsed / 3600))
minutes=$(( (seconds_elapsed % 3600) / 60 ))
seconds=$((seconds_elapsed % 60))
echo "Time elapsed: ${hours}h:${minutes}m:${seconds}s"