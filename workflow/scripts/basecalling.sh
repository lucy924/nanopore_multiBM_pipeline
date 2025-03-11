#!/bin/bash
set -e

# snakemake rule
# rule basecall_sample:
#     # To one bam file hopefully
#     input:
#         sample_pod5_dir = f"config/{SAMPLE}/"
#     output:
#         sample_raw_bam = f"results/basecalled/{SAMPLE}.unaligned.bam"
#     params:
#         sample_name = f"{SAMPLE}",
#         nf_threads = config['nextflow']['threads'],
#         nf_profile = config['nextflow']['profile']
#     log: 
#         "logs/basecall_sample.log",
#     script:
#         "scripts/basecalling.sh"

# Never actually got this script working

POD5_DIR=${snakemake_input[sample_pod5_dir]}
BAM_OUTPUT=${snakemake_input[sample_raw_bam]}
SAMPLE_NAME=${snakemake_params[sample_name]}
THREADS=${snakemake_params[nf_threads]}
PROFILE=${snakemake_params[nf_profile]}
LOG=${snakemake_log}

echo "POD5 directory: ${POD5_DIR}" >> $LOG
echo "BAM output file: ${BAM_OUTPUT}" >> $LOG
echo "Sample name: ${SAMPLE_NAME}" >> $LOG

# sample_name="LPFT"
# input_pod5s="/external/data/lucy/20241107-BiomarkerProject/late_pod5s"
# out_dir="${wd}/late_bams"

out_dir="$(dirname "${BAM_OUTPUT}")" 

basecaller_cfg="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
remora_cfg="dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v2.0.1"

nextflow run epi2me-labs/wf-basecalling \
    --sample_name "${SAMPLE_NAME}.unaligned" \
    --input "${POD5_DIR}" \
    --out_dir "${out_dir}" \
    --output_fmt bam \
    --basecaller_cfg "${basecaller_cfg}" \
    --remora_cfg "${remora_cfg}" \
    --basecaller_basemod_threads "${THREADS}" \
    -profile "${PROFILE}"

# ---------------------------------------------------------- #
# Additional basecaller info in case my stuff doesn't work

# - dna_r10.4.1_e8.2_260bps_sup@v4.1.0_5mCG_5hmCG@v2
# - dna_r10.4.1_e8.2_400bps_sup@v4.1.0_5mCG_5hmCG@v2
# - dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1
# - dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mC_5hmC@v1
# - dna_r10.4.1_e8.2_400bps_sup@v4.3.0_6mA@v2
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v2
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v2.0.1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1
# - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v2
