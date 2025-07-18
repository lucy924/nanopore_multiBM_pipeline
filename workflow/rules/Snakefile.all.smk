# Making the dag:  
# snakemake --dag [your target file] | dot -Tsvg > dag.svg
# snakemake --dag make_minknow_targets | dot -Tsvg > dag.snv_run.svg
# snakemake --dag small_data_test | dot -Tsvg > dag.small_data_test.svg
# snakemake --dag run_methylCS | dot -Tsvg > dag.run_methylCS.svg
# snakemake --dag test_things | dot -Tsvg > dag.test_things.svg
# snakemake test_things --use-conda

# import snakemake

report: "reports/workflow.rst"
configfile: "config/config.yaml"

include: "panel_prep.smk"
# include: "mod_prep.smk"

# SAMPLE = config['sample']

cores = config.get("max_cores", 16)  # Default to 16 if not set
# snakemake.snakemake("Snakefile", cores=cores)

rule test_things:
    # snakemake --dag test_things | dot -Tsvg > dag.test_things.svg
    # snakemake test_things --use-conda
    input:
        panel_results = f"results/{SAMPLE}/{SAMPLE}.panel_results.csv"

rule run_sample_panel_for_BM_classifier:
    # snakemake --dag run_sample_panel_for_BM_classifier | dot -Tsvg > dag.run_sample_panel_for_BM_classifier.svg
    # snakemake run_sample_panel_for_BM_classifier --use-conda
    input:
        panel_results = f"results/{SAMPLE}/{SAMPLE}.panel_results.csv"

rule get_input_for_minknow:
    input:
        "minknow_input/targets.minknow.bed"

rule merge_bams:
    input:
        config['bam_pass_directory']
    output:
        # f"results/{SAMPLE}/{SAMPLE}.bam"
        os.path.abspath(f"results/{SAMPLE}/{SAMPLE}.bam")
    params:
        path_temp = os.path.abspath(f"results/{SAMPLE}/temp"),
        threads = cores
    log:
        f"logs/{SAMPLE}/merge_bams.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.merge_bams.tsv"
    script:
        "scripts/merge_bams.sh"
        # "samtools merge -o {output.sample_bam} {input.bam_pass_dir}/*.bam 2> {log}"

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
        "scripts/run_epi2me_humvar.sh"

rule bgzip_clinvar_vcf:
    input:
        vcf_clinvar = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf"
    output:
        vcf_clinvar_gz = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.bgzip_clinvar_vcf.tsv"
    shell:
        "bgzip {input}"

rule tabix_clinvar_vcf_gz:
    input:
        vcf_clinvar_gz = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz"
    output:
        vcf_clinvar_gz_tbi = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz.tbi"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.tabix_clinvar_vcf_gz.tsv"
    shell:
        "tabix {input}"

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
        "scripts/snv_annotation.py"

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
        "scripts/sv_annotation.py"


rule combine_bedmethyls:
    input:
        bedmeth1 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.1.bedmethyl.gz",
        bedmeth2 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.2.bedmethyl.gz",
        bedmeth3 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.ungrouped.bedmethyl.gz"
    output:
        temp(f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed")
    params:
        threads = config['threads']
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.combine_bedmethyls.tsv"
    shell:
        (
        "bgzip -dc {input.bedmeth1} {input.bedmeth2} {input.bedmeth3} | "
        "sort --parallel={params.threads} -k1,1 -k2,2n > {output}"
        )

rule convert_bedmethyl_to_DSS:
    input: 
        f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed"
    output: 
        f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.convert_bedmethyl_to_DSS.tsv"
    script:
        "scripts/convert_DSS.py"


# rule transform_DSS_for_ImIn:
    #     input:
    #         dss_file = f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    #     params:
    #         sample_name = SAMPLE

        # """
        # awk -v OFS='\t' 'BEGIN{print "chr","pos","N","X"}{print $1,$2,($12+$13),$13}' {input} > {output}
        # """


rule prep_for_getting_betas:
    input:
        IlluminaEPIC_genomic_locations_hg38 = "resources/IlluminaEPIC_genomic_locations_hg38.csv",
        dss_file = f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    output:
        pre_beta = f"results/{SAMPLE}/mod_calling/{SAMPLE}.pre_beta.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.prep_for_betas.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.prep_for_getting_betas.tsv"
    script:
        "scripts/process_result_before_betas.py"

rule add_betas:
    input:
        pre_beta = f"results/{SAMPLE}/mod_calling/{SAMPLE}.pre_beta.csv"
    output:
        # post_beta = temp(f"results/{SAMPLE}/mod_calling/{SAMPLE}.post_beta.csv")
        post_beta = f"results/{SAMPLE}/mod_calling/{SAMPLE}.post_beta.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.add_betas.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.add_betas.tsv"
    conda: 
        "envs/methylcibersort.yaml"
    script:
        "scripts/add_betas.R"

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
        "scripts/modification_calling.py"

rule run_methylCS:
    # script looks a bit weird because the program does things with the names to make the final files
    input:
        beta_matrix = f"results/{SAMPLE}/mod_calling/{SAMPLE}.methatlas.csv"
    output:
        cibersort_mix_matrix_upload_file = f"results/{SAMPLE}/methylCS/{SAMPLE}.CS_mix_matrix.txt",
        cibersort_bladder_ref_upload_file = f"results/{SAMPLE}/methylCS/{SAMPLE}.CS_bladder_ref.txt"
    params:
        cibersort_mix_matrix_upload_fname = f"results/{SAMPLE}/methylCS/{SAMPLE}.CS_mix_matrix",
        sample_name = {SAMPLE}
    conda: 
        "envs/methylcibersort.yaml"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.run_methylCS.tsv"
    shell:
        """
        if ! Rscript -e "if (!requireNamespace('MethylCIBERSORT', quietly=TRUE)) quit(status=1)"; then
            R CMD INSTALL workflow/scripts/MethylCIBERSORT/v2.0/MethylCIBERSORT_0.2.0.tar.gz
        fi

        Rscript workflow/scripts/methylcibersort.r {input.beta_matrix} {params.cibersort_mix_matrix_upload_fname} {output.cibersort_bladder_ref_upload_file} {params.sample_name}
        """

rule run_CIBERSORTX:
    input:
        cibersort_mix_matrix_upload_file = f"results/{SAMPLE}/methylCS/{SAMPLE}.CS_mix_matrix.txt",
        cibersort_bladder_ref_upload_file = f"results/{SAMPLE}/methylCS/{SAMPLE}.CS_bladder_ref.txt",
        standin = "resources/test3.CS_bladder_ref.csv"
    output:
        cibersortx_output = f"results/{SAMPLE}/methylCS/CIBERSORTx_{SAMPLE}_Results.csv"
    params:
        username = config["cibersortx"]["username"],
        token = config["cibersortx"]["token"]
    shell:
        ("cp {input.standin} {output.cibersortx_output}")

# rule run_methatlas:
    #     # find immune infiltrate of tumour
    #     # TCGA data currently being used
    #     input:
    #         data = f"results/{SAMPLE}/mod_calling/{SAMPLE}.methatlas.csv",
    #         bladder_atlas = "resources/test3.CS_bladder_ref.csv"
    #     output:
    #         f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.methatlas_deconv_output.csv"
    #     params:
    #         out_dir = f"results/{SAMPLE}/immune_infiltrate"
    #     conda: 
    #         "envs/methatlas.yaml"
    #     benchmark:
    #         f"logs/{SAMPLE}/{SAMPLE}.benchmark.run_methatlas.tsv"
    #     shell:
    #         "python meth_atlas/deconvolve.py -a {input.bladder_atlas} --out_dir {params.out_dir} {input.data}"


# rule immune_infiltrate_methatlas:
    #     input:
    #         immune_reference_dataset="resources/ref_atlas_bladder.csv",
    #         epic_probe_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.methatlas.csv",
    #         methatlas_deconvolution = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.methatlas_deconv_output.csv"
    #     output:
    #         results_json_imin = f"results/{SAMPLE}/results.imin.json",
    #         immune_results = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_infiltrate.csv"
    #     log:
    #         f"logs/{SAMPLE}/{SAMPLE}.immune_infiltrate.log"
    #     script:
    #         "scripts/get_immune_infiltrate.py"


rule immune_infiltrate_mCS:
    input:
        panel_metadata = "config/panel_metadata.csv",
        mCS_results = f"results/{SAMPLE}/methylCS/CIBERSORTx_{SAMPLE}_Results.csv"
    output:
        immune_results = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.immune_infiltrate.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.immune_infiltrate_mCS.tsv"
    script:
        "scripts/get_immune_infiltrate.mCS.py"

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
        "scripts/collate_results_for_preclintrial.py"

rule get_patient_results_csv:
    input:
        panel_metadata = "config/panel_metadata.csv",
        snv_results = f"results/{SAMPLE}/snv_annotation/{SAMPLE}.snv_results.csv",
        sv_results = f"results/{SAMPLE}/sv_annotation/{SAMPLE}.sv_results.csv",
        mod_results = f"results/{SAMPLE}/mod_calling/{SAMPLE}.methatlas.csv",
        immune_results = f"results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
    output:
        patient_results = f"results/{SAMPLE}/{SAMPLE}.results.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.get_patient_results.log"
    benchmark:
        f"logs/{SAMPLE}/{SAMPLE}.benchmark.get_patient_results.tsv"
    script:
        "scripts/get_patient_results.py"


rule get_scores:
    input:
        panel_metadata = "config/panel_metadata.scores.csv",
        combined_panel_results = f"results/{SAMPLE}/{SAMPLE}.all_panel_results.csv"
    output:
        scores = f"results/{SAMPLE}/{SAMPLE}.scores.csv"
    log:
        f"logs/{SAMPLE}/{SAMPLE}.get_scores.log"
    script:
        "scripts/get_scores.py"
        

rule generate_report:
    # find immune infiltrate of tumour
    #TODO: take out ratios from the imin.json input
    input:
        methatlas_deconv_fig = f"{SAMPLE}.methatlas_deconv_plot.png",
        results_json_imin = f"results/{SAMPLE}/results.imin.json"
    output:
        "file"
    script:
        "scripts/generate_report_template.py"
