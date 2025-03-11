rule index_ref:
    input:
        "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
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
        immune_reference_dataset="resources/ref_atlas_bladder.csv",
        epic_locs_hg38="resources/IlluminaEPIC_genomic_locations_hg38.csv",
    output:  # TODO: something has gone weird after I added the genebuffed section, it's not getting enough coverage anymore. suspect it's passing the wrong version of the panel somewhere.
        # panel_genebuffed_csv = "results/panel_genebuffed.csv",  # this has promoter+downstream regions for whole gene regions
        panel_bed="results/panel.bed",  # Note this has targeted buffer regions added
        all_targets="results/panel_prep/all_targets.bed",  # Note this has targeted buffer regions added and immune infiltrate methylation locations
        panel_csv_b4_targets="debug/panel.csv",
        panel_csv_targets="debug/panel_targets.csv"
    log:
        "logs/make_minknow_input/make_panel_bed.log",
    script:
        "scripts/make_panel_bed.py"

rule run_make_adaptive_ref:
    input:
        all_targets="results/panel_prep/all_targets.bed",
        ref_fasta="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        fasta_index="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        chrom_sizes=f"resources/hg38_no_alt.chrom_sizes",
    output:
        # temp(f"results/output_{config['flanking_bp']}.txt"),
        minknow_target_bed=f"minknow_input/targets.minknow.{config['flanking_bp']}.bed", # This has bigger buffer regions for minknow input
        minknow_target_fasta=f"minknow_input/targets.minknow.{config['flanking_bp']}.fasta",
        sorted_targets=f"results/make_minknow_input/sorted_all_targets.{config['flanking_bp']}.bed",
        ini_targets=f"results/make_minknow_input/ini_all_targets.{config['flanking_bp']}.bed"
    log:
        f"logs/make_minknow_input/all.{config['flanking_bp']}.log"
    script:
        "scripts/make_adaptive_ref.sh"

rule check_coverage:
    input:
        # f"results/output_{config['flanking_bp']}.txt",
        minknow_target_bed=f"minknow_input/targets.minknow.{config['flanking_bp']}.bed",
    output:
        # "results/final_output.txt",
        final_bed_name = "minknow_input/targets.minknow.bed"  # This has bigger buffer regions for minknow input
    params:
        min_cov = config["min_genome_coverage"],
        max_cov = config["max_genome_coverage"]
    log:
        "logs/make_minknow_input/check_criteria.log"
    script:
        "scripts/calculate_coverage.py"
