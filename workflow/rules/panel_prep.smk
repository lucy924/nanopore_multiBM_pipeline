PROJECT = config['project_name']
BUFFER = config['buffersize_bp']

rule index_ref:
    input:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",  # TODO: make this an input in the config file
    output:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
    shell:
        "samtools faidx {input}"

rule get_chrom_sizes:
    input:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
    output:
        "resources/hg38_no_alt.chrom_sizes"
    shell:
        "cut -f1,2 {input} > {output}"

rule make_panel_bed:
    input:
        panel_csv=config["panel_metadata"],
        # immune_reference_dataset="resources/test3.CS_bladder_ref.csv",  # this selection is too small for adaptive seq coverage
        immune_reference_dataset="resources/ref_atlas_bladder.csv",
        epic_locs_hg38="resources/IlluminaEPIC_genomic_locations_hg38.csv",
    output:
        panel_bed=f"results/{PROJECT}/minknow_input_supp/biomarker_panel.bed",
        all_targets=f"results/{PROJECT}/minknow_input/targets.bed",
    log:
        "logs/make_minknow_input/make_panel_bed.log",
    script:
        "../scripts/make_panel_bed.py"

rule run_make_adaptive_ref:
    input:
        all_targets=f"results/{PROJECT}/minknow_input/targets.bed",
        ref_fasta="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        fasta_index="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        chrom_sizes=f"resources/hg38_no_alt.chrom_sizes",
    output:
        # temp(f"results/output_{config['flanking_bp']}.txt"),
        minknow_target_bed=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{BUFFER}.bed",
        minknow_target_fasta=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{BUFFER}.fasta",
        sorted_targets=f"results/{PROJECT}/minknow_input_supp/sorted_all_targets.{BUFFER}.bed",
        ini_targets=f"results/{PROJECT}/minknow_input_supp/ini_all_targets.{BUFFER}.bed"
    log:
        f"logs/{PROJECT}/minknow_input_supp/make_adaptive_ref.{BUFFER}.log"
    script:
        "../scripts/make_adaptive_ref.sh"

rule check_coverage:
    input:
        minknow_target_bed=f"results/{PROJECT}/minknow_input_supp/targets.minknow.{BUFFER}.bed"
    output:
        final_bed_name = f"results/{PROJECT}/minknow_input/targets.buffed.bed"
    params:
        min_cov = config["min_genome_coverage"],
        max_cov = config["max_genome_coverage"],
        buffer = BUFFER
    log:
        f"logs/{PROJECT}/minknow_input_supp/check_criteria.{BUFFER}.log"  # For some reason this errors
    script:
        "../scripts/calculate_coverage.py"
