#!/bin/bash
#SBATCH --job-name=3_phase
#SBATCH --output=3_phase-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=30G
#SBATCH --cpus-per-task=40
#SBATCH --array=22
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

## SHAPEIT5_v2 changes:
# using only ['AMR', 'CAH', 'EUR', 'AAC', 'AFR'] samples (/cluster/home/seger/scripts/redlat_tangl/global_ancestry/genotools_include_anc.py)
# using only trios & duos that were confirmed by KING/PC-Relate (on kosik)

ml load bcftools
ml load cluster/singularity/3.11.0

# get threads
threads=${SLURM_CPUS_PER_TASK}

# get chromosome
CHR=${SLURM_ARRAY_TASK_ID}
echo "chr${CHR}"

## directories
MAP_DIR="/cluster/home/seger/software/shapeit5/resources/maps/b38"
PED_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/GENESIS_v2"
INDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5_v3/2_sample_qc"
OUTDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/SHAPEIT5_v3/3_phase_common"
mkdir -p ${OUTDIR}

## files
MAP="${MAP_DIR}/chr${CHR}.b38.gmap.gz"
PED="${PED_DIR}/redlat_tangl_confirmed.ped"
IN="${INDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.include_anc.bcf"
OUT="${OUTDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.include_anc.maf0001.shapeit5_common.bcf"
LOG="${OUTDIR}/redlat_tangl.chr${CHR}.qceed.snps.mend.include_anc.maf0001.shapeit5_common.log"

start=$(date +%s)

# phase_common
singularity run \
--bind /usr/lib/locale/ \
--bind /cluster/:/cluster/ \
/cluster/home/seger/software/shapeit5_v5.1.1.sif \
phase_common_static \
    --filter-maf 0.001 \
    --region chr${CHR} \
    --input ${IN} \
    --pedigree ${PED} \
    --map ${MAP} \
    --output ${OUT} \
    --thread ${threads} \
    --log ${LOG}

bcftools index -f ${OUT}

end=$(date +%s)

## print time elapsed in hours, minutes, and seconds
seconds_elapsed=$((end - start))
hours=$((seconds_elapsed / 3600))
minutes=$(( (seconds_elapsed % 3600) / 60 ))
seconds=$((seconds_elapsed % 60))
echo "Time elapsed: ${hours}h:${minutes}m:${seconds}s"