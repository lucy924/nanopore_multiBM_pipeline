#!/bin/bash -e
# script name: make_adaptive_ref.sh
# prepare reference file for adaptive sampling on ONT

# ----------------------- #
# Original script by
# Stephane Plaisance (VIB-NC) 2021/03/23; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

# Modified by Lucy Picard for biomarker pipeline
# ----------------------- #

version="1.1, 2025_02_08 - LP"

REF=${snakemake_input[ref_fasta]}
BED=${snakemake_input[all_targets]}
CHROM_SIZES=${snakemake_input[chrom_sizes]}
opts=${snakemake_params[buffersize_bp]}
LOG=${snakemake_log}

# check executable present
declare -a arr=( "samtools" "bedtools" "gawk")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# additional side length
BASES_TO_EXPAND_PER_SIDE=${opts:-5000}

# OUTPUTS
SLOPPED_BED=${snakemake_output[minknow_target_bed]}
SUBSETTED_FASTA=${snakemake_output[minknow_target_fasta]}

echo "BED file: ${BED}" >> $LOG
echo "Output file: ${snakemake_output[sorted_targets]}" >> $LOG
echo "sort input BED by chr then start" >> $LOG
sort -k 1V,1 -k 2n,2 ${BED} \
  > ${snakemake_output[sorted_targets]}

echo "create expanded BED" >> $LOG
 
echo ${BASES_TO_EXPAND_PER_SIDE}
echo ${BASES_TO_EXPAND_PER_SIDE}
echo ${snakemake_output[sorted_targets]}
echo ${CHROM_SIZES}

bedtools slop -l 2000 -r 2000 -i /home/dejlu879/ProjectProtocol/nanopore_multiBM_pipeline/results/BCG_on_NMIBC/minknow_input_supp/sorted_all_targets.2000.bed -g /home/dejlu879/ProjectProtocol/nanopore_multiBM_pipeline/resources/hg38_no_alt.chrom_sizes > /home/dejlu879/ProjectProtocol/nanopore_multiBM_pipeline/results/BCG_on_NMIBC/minknow_input_supp/ini_all_targets.test.bed

bedtools slop \
  -l ${BASES_TO_EXPAND_PER_SIDE} \
  -r ${BASES_TO_EXPAND_PER_SIDE} \
  -i ${snakemake_output[sorted_targets]} \
  -g ${CHROM_SIZES} \
  > ${snakemake_output[ini_targets]}

# echo "merge region overlaps where present and collapse their descriptions as a csv list" >> $LOG
# bedtools merge -i ${snakemake_output[ini_targets]} \
#   -c 4 \
#   -o collapse \
#   > ${SLOPPED_BED}

# echo "print total reference width" >> $LOG
# TOT_WIDTH=$(gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' \
#   ${SLOPPED_BED})
# echo "# total reference width in ${SLOPPED_BED} is $TOT_WIDTH bps" >> $LOG
# if [ $TOT_WIDTH < 500 ]; then
#   exit 1
# fi

# echo "extract fasta sequences" >> $LOG
# bedtools getfasta -fi ${REF} \
#   -bed ${SLOPPED_BED} \
#   -fo ${SUBSETTED_FASTA} \
#   -name


# # This is the final file which you will upload into MinKNOW:
# echo "# the file ${SUBSETTED_FASTA} can be used in Minknow for adaptive sequencing" >> $LOG
# echo "# the file ${SUBSETTED_FASTA} can be used in Minknow for adaptive sequencing"

# exit 0
