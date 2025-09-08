#!/bin/bash

# Variant positions
VAR1=44908684
VAR2=44908822

# directories
# DIR1="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/LAI/RFmix/Tractor"
DIR1="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/LAI/RFmix_v2/Tractor"
DIR2="/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/LAI/APOE"
mkdir -p ${DIR2}
cd ${DIR2}

# File list
# FILES=(
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc0.dosage.txt.gz
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc0.hapcount.txt.gz
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc1.dosage.txt.gz
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc1.hapcount.txt.gz
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc2.dosage.txt.gz
#   ${DIR1}/redlat_tangl.AFR_AMR_EUR.chr19.qceed.miss005.maf0001.shapeit5_common.anc2.hapcount.txt.gz
# )

FILES=(
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc0.dosage.txt.gz
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc0.hapcount.txt.gz
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc1.dosage.txt.gz
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc1.hapcount.txt.gz
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc2.dosage.txt.gz
  ${DIR1}/redlat_tangl.chr19.qceed.include_anc.maf0001.shapeit5_common.anc2.hapcount.txt.gz
)

for FILE in "${FILES[@]}"; do
  echo "Processing $FILE..."

  # Extract prefix for output file naming
  ANC=$(echo "$FILE" | grep -o 'anc[0-9]')
  TYPE=$(echo "$FILE" | grep -oE '(dosage|hapcount)')

  # OUTFILE="redlat_tangl.AFR_AMR_EUR.APOE.${ANC}.${TYPE}.txt"
  OUTFILE="redlat_tangl.include_anc.APOE.${ANC}.${TYPE}.txt"

  # Extract lines: header and the two desired variants
  zcat "$FILE" | awk -v var1="$VAR1" -v var2="$VAR2" 'NR==1 || $2==var1 || $2==var2' > tmp_selected.txt

  # Extract variant names for new column headers
  header1=$(awk -v var1="$VAR1" '$2==var1 {print $1":"$2":"$4":"$5}' tmp_selected.txt)
  header2=$(awk -v var2="$VAR2" '$2==var2 {print $1":"$2":"$4":"$5}' tmp_selected.txt)

  # Remove first 5 columns and transpose the remaining
  cut -f6- tmp_selected.txt | \
  awk '
  {
    for (i=1; i<=NF; i++) {
      a[NR,i] = $i
    }
    max_i = (NF > max_i ? NF : max_i)
    max_NR = NR
  }
  END {
    for (i=1; i<=NF; i++) {
      for (j=1; j<=NR; j++) {
        printf "%s%s", a[j,i], (j<NR ? OFS : ORS)
      }
    }
  }
  ' OFS='\t' > tmp_transposed.txt

  # Add header
  ( echo -e "SAMPLE\t$header1\t$header2"; cat tmp_transposed.txt ) > "$OUTFILE"

  echo "Wrote output to $OUTFILE"
done

# Clean up
rm -f tmp_selected.txt tmp_transposed.txt
