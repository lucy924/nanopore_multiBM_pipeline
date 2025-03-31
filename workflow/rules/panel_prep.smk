PROJECT = config['project_name']

rule index_ref:
    input:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
    # conda: 
    #     "../envs/general.yaml"
    shell:
        "samtools faidx {input}"

rule get_chrom_sizes:
    input:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
    output:
        "resources/hg38_no_alt.chrom_sizes"
    # conda: 
    #     "../envs/general.yaml"
    shell:
        "cut -f1,2 {input} > {output}"

rule make_panel_bed:
    input:
        panel_csv=config["panel_metadata"],
        immune_reference_dataset="resources/ref_atlas_bladder.csv",
        epic_locs_hg38="resources/IlluminaEPIC_genomic_locations_hg38.csv",
    output:
        panel_bed=f"results/{PROJECT}/minknow_input_supp/biomarker_panel.bed",
        all_targets=f"results/{PROJECT}/minknow_input/targets.bed",
    log:
        f"logs/{PROJECT}/minknow_input_supp/make_panel_bed.log",
    # conda: 
    #     "../envs/general.yaml"
    script:
        "../scripts/make_panel_bed.py"

rule run_make_adaptive_ref:
    input:
        all_targets=f"results/{PROJECT}/minknow_input/targets.bed",
        ref_fasta="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        fasta_index="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        chrom_sizes=f"resources/hg38_no_alt.chrom_sizes",
    output:
        # temp(f"results/output_{config['buffersize_bp']}.txt"),
        minknow_target_bed=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{config['buffersize_bp']}.bed",
        minknow_target_fasta=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{config['buffersize_bp']}.fasta",
        sorted_targets=f"results/{PROJECT}/minknow_input_supp/sorted_all_targets.{config['buffersize_bp']}.bed",
        ini_targets=f"results/{PROJECT}/minknow_input_supp/ini_all_targets.{config['buffersize_bp']}.bed"
    log:
        f"logs/{PROJECT}/minknow_input_supp/make_adaptive_ref.{config['buffersize_bp']}.log"
    # conda: 
    #     "../envs/general.yaml"
    script:
        "../scripts/make_adaptive_ref.sh"

rule check_coverage:
    input:
        # f"results/output_{config['buffersize_bp']}.txt",
        minknow_target_bed=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{config['buffersize_bp']}.bed",
    output:
        # "results/final_output.txt",
        final_bed_name = f"results/{PROJECT}/minknow_input/targets.buffed.bed"
    params:
        min_cov = config["min_genome_coverage"],
        max_cov = config["max_genome_coverage"]
    log:
        f"logs/{PROJECT}/minknow_input_supp/check_criteria.log"
    # conda: 
    #     "../envs/general.yaml"
    script:
        "../scripts/calculate_coverage.py"
