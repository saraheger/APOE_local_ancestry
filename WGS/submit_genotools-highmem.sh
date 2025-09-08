#!/bin/bash
#SBATCH --job-name=genotools
#SBATCH --output=genotools-%A.out
#SBATCH --time=48:00:00
#SBATCH --partition=highmem
#SBATCH --mem=600G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eger@ucsb.edu

source ~/.bashrc
ml load cluster/plink/1.90
ml load cluster/singularity

conda activate /cluster/home/jtaylor/micromamba/envs/genotools_wgs

## directories
MAINDIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25"
INDIR="${MAINDIR}/global_ancestry"
OUTDIR="${INDIR}/genotools"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

# this doesn't exclude the miss005 samples - plink file was incorrectly named
genotools \
    --singularity \
    --bfile ${INDIR}/redlat_tangl.auto.qceed.miss005 \
    --out ${OUTDIR}/redlat_tangl.genotools \
    --ancestry \
    --ref_panel /cluster/home/seger/references/genotools/ref_panel/1kg_30x_hgdp_ashk_ref_panel \
    --ref_labels /cluster/home/seger/references/genotools/ref_panel/1kg_30x_hgdp_ashk_ref_panel_labels.txt