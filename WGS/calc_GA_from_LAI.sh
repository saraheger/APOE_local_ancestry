#!/bin/bash

# srun -p interactive --pty --time=24:00:00 --ntasks 2 bash

ml load python

ADMIXTURE_DIR="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/global_ancestry/ADMIXTURE"
DIR1="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/LAI/RFmix_v2"
DIR2="${DIR1}/PlotCoversion"
DIR3="${DIR1}/CollapseAncestry"
DIR4="${DIR1}/EstimateGlobalAncestry"

id_list="${ADMIXTURE_DIR}/include_anc_sampleIDs.txt"

mkdir -p ${DIR2}
cd ${DIR2}

### Plot Conversion
# code to convert RFmixv2 msp output to a format to plot up painted karyograms.
# https://github.com/eatkinson/ancestry_pipeline/blob/
# 218b6da4c5f6155d0d40af597ef9ce8c726aa78e/RFmixv2_plotconversion.txt

# First make a file with snp locations
for i in {1..22}; do \
    PREFIX="redlat_tangl.chr${i}.qceed.include_anc.maf0001.shapeit5_common.LAI"
    awk '{ if ($1!~"#") print $0 }' ${DIR1}/${PREFIX}.msp.tsv | \
    awk '{ if ($1!~"#") print $0 }'| \
    cut -f2,4 > ${DIR2}/${PREFIX}.snp_loc ; \
    done
 
# Process msp output also to match original RFmix Viterbi format 
for i in {1..22}; do \
    PREFIX="redlat_tangl.chr${i}.qceed.include_anc.maf0001.shapeit5_common.LAI"
    awk '{ if ($1!~"#") print $0 }' ${DIR1}/${PREFIX}.msp.tsv | \
    awk '{ if ($1!~"#") print $0 }' > ${DIR2}/${PREFIX}.msp1.tsv ; \
    done

for i in {1..22}; do \
    PREFIX="redlat_tangl.chr${i}.qceed.include_anc.maf0001.shapeit5_common.LAI"
    cut -f7- ${DIR2}/${PREFIX}.msp1.tsv > ${DIR2}/${PREFIX}.Viterbi.tsv ; \
    done

### Collapse Ancestry
# Now you can run Alicia's collapsing ancestry script 
# as normal to create bed files for each individual, 
# putting in a list with the indiv IDs 

mkdir -p ${DIR3}
cd ${DIR3}

# does all the chromosome at once
awk '{print $1}' ${id_list} | while read id; do
    outfile=${DIR3}/${id}.unkinked.iterated.B.bed
    if [ ! -f "$outfile" ]; then
        python /cluster/home/seger/software/ancestry_pipeline_eatkinson/collapse_ancestry1.py \
            --rfmix ${DIR2}/redlat_tangl.chr1.qceed.include_anc.maf0001.shapeit5_common.LAI.Viterbi.tsv \
            --snp_locations ${DIR2}/redlat_tangl.chr1.qceed.include_anc.maf0001.shapeit5_common.LAI.snp_loc \
            --ind $id \
            --ind_info ${id_list} \
            --pop_labels "AFR,AMR,EUR" \
            --out ${DIR3}/$id.unkinked.iterated
    fi
done

### Estimate Global Ancestry
mkdir -p ${DIR4}
cd ${DIR4}

# create bed_list
awk -v d3="${DIR3}" '{print d3"/"$1".unkinked.iterated.A.bed\t"d3"/"$1".unkinked.iterated.B.bed"}' \
    ${id_list} > ${DIR4}/redlat_tangl.bed_list.txt

# Estimate global ancestry proportions from local ancestry inference
python /cluster/home/seger/software/ancestry_pipeline_eatkinson/lai_global.py \
    --bed_list ${DIR4}/redlat_tangl.bed_list.txt \
    --ind_list ${id_list} \
    --pops 'AFR,AMR,EUR' \
    --out ${DIR4}/redlat_tangl.lai_global.txt

