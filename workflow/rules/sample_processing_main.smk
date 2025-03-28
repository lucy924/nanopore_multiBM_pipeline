
SAMPLE = config['sample']
cores = config.get("max_cores", 16)  # Default to 16 if not set


rule run_wf_humvar:
    input:
        # sample_bam = os.path.abspath(f"results/{SAMPLE}/{SAMPLE}.bam"),
        sample_bam = config['bam_pass_directory'],
        reference = os.path.abspath("resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"),
        targets_bed_file = os.path.abspath("results/panel.bed"),
        tandem_repeat_bed = os.path.abspath("resources/hg38.trf.bed.gz")
    output:  #TODO: add all necessary outputs
        vcf_clinvar = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf",
        vcf_all = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp.vcf.gz",
        vcf_sv_gz = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_sv.vcf.gz",
        mod_strand1 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.1.bedmethyl.gz",
        mod_strand2 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.2.bedmethyl.gz",
        mod_unstranded = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.ungrouped.bedmethyl.gz",
        panel_coverage_data = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.gene_summary.tsv",
        flag = os.path.abspath(f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf-humvar.finished.flag")
    params:
        sample_name = SAMPLE,
        bam_min_coverage = 0,  # will be 20 in final
        nf_basecaller = config["nextflow"]['basecaller'],
        nf_threads = config["nextflow"]['threads'],
        nf_profile = config["nextflow"]['profile'],
        abs_path_to_results = os.path.abspath("results")
    log:
        f"logs/{SAMPLE}/run_wf_humvar.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.run_wf_humvar.tsv"
    script:
        "../scripts/run_epi2me_humvar.sh"
        

rule snv_annotation:
    # SNV annotation by ClinVar
    input:
        panel_metadata = "config/panel_metadata.csv",
        vcf_clinvar_gz = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz",
        vcf_clinvar_gz_tbi = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz.tbi",
        vcf_all = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp.vcf.gz",
    output:
        snv_csv = f"results/{SAMPLE}/snv_annotation/{SAMPLE}.raw_snv_results.csv",
        snv_panel_csv = f"results/{SAMPLE}/snv_annotation/{SAMPLE}.snv_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.snv_annotation.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.snv_annotation.tsv"
    script:
        "../scripts/snv_annotation.py"


rule sv_annotation:
    # Structural variant (SV) calling by SNPEff – identifies repeat expansions
    # TODO: have not yet figured out svs, have none in my current panel. Do have some that were sequenced though, so can make stuff up to test in the future.
    input:
        panel_metadata = "config/panel_metadata.csv",
        vcf_sv_gz = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_sv.vcf.gz",
    output:
        sv_csv = f"results/{SAMPLE}/sv_annotation/{SAMPLE}.raw_sv_results.csv",
        sv_panel_csv = f"results/{SAMPLE}/sv_annotation/{SAMPLE}.sv_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.sv_annotation.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.sv_annotation.tsv"
    script:
        "../scripts/sv_annotation.py"


rule modification_calling:
    # Methylation identification - will take at least 15 mins
    #TODO: can probs make it faster by combining get_results_df with downstream processes
    input:
        panel_metadata = "config/panel_metadata.csv",
        post_beta = f"results/{SAMPLE}/mod_calling/{SAMPLE}.post_beta.csv"
    output:
        epic_probe_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.methatlas.csv",
        panel_mod_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.mod_results.csv",
        panel_rawmod_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.rawmod_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.modification_calling.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.modification_calling.tsv"
    script:
        "../scripts/modification_calling.py"

rule immune_infiltrate_mCS:
    input:
        panel_metadata = "config/panel_metadata.csv",
        mCS_results = f"results/{SAMPLE}/methylCS/CIBERSORTx_{SAMPLE}_Results.txt"
    output:
        immune_results = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.immune_infiltrate.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.immune_infiltrate_mCS.tsv"
    script:
        "../scripts/get_immune_infiltrate.mCS.py"


rule collate_results_for_preclin_trial:
    input:
        panel_metadata = "config/panel_metadata.csv",
        snv_results = f"results/{SAMPLE}/snv_annotation/{SAMPLE}.snv_results.csv",
        mod_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.mod_results.csv",
        immune_results = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
    output:
        panel_results = f"results/{SAMPLE}/{SAMPLE}.panel_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.collate_results_for_preclin_trial.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.collate_results_for_preclin_trial.tsv"
    script:
        "../scripts/collate_results_for_preclintrial.py"


