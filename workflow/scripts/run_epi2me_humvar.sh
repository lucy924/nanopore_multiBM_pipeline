#!/bin/bash -e
# script name: run_epi2me_humvar.sh

SAMPLE_BAM=${snakemake_input[sample_bam]}
REFERENCE=${snakemake_input[reference]}
TARGETS=${snakemake_input[targets_bed_file]}
TANDEM_REPEAT=${snakemake_input[tandem_repeat_bed]}
SAMPLE=${snakemake_params[sample_name]}
BAM_MIN_COV=${snakemake_params[bam_min_coverage]}
THREADS=${snakemake_params[nf_threads]}
PROFILE=${snakemake_params[nf_profile]}
RESULTS_PATH=${snakemake_params[abs_path_to_results]}
LOG=${snakemake_log}

BASECALLER=${snakemake_params[nf_basecaller]}  # shouldn't need this in final version when everything is run at once and I don't need to combine different runs to get a good readc depth 

cd results/${SAMPLE}/wf-humvar/
# mkdir -p results
# cp "${TARGETS}" results/.

nextflow run epi2me-labs/wf-human-variation \
    --bam "${SAMPLE_BAM}" \
    --ref "${REFERENCE}" \
    --bed "${TARGETS}" \
    --out_dir . \
    --sample_name "${SAMPLE}" \
    --sv \
    --snp \
    --mod \
    --str \
    --phased \
    --override_basecaller_cfg "${BASECALLER}"\
    --output_gene_summary \
    --output_xam_fmt bam \
    --modkit_args "--preset traditional" \
    --bam_min_coverage "${BAM_MIN_COV}" \
    --threads "${THREADS}" \
    --ubam_map_threads "${THREADS}" \
    --ubam_sort_threads "${THREADS}" \
    -profile "${PROFILE}"

# echo "removing file:"
# echo "results/$(basename ${TARGETS})"
# rm "results/$(basename ${TARGETS})"

# Testing moving files
# echo "mv results/test_file_in_wf-humvar.txt ${RESULTS_PATH}/."
# touch results/test_file_in_wf-humvar.txt
# mv results/test_file_in_wf-humvar.txt ${RESULTS_PATH}/.
# mv results/*.bedmethyl.gz ${RESULTS_PATH}/.

touch ${snakemake_output[flag]}

# Not using below parameters for now
# --tr_bed {input.tandem_repeat_bed} \
# --cnv \
# --use_qdnaseq \
