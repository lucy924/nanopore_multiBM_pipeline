
SAMPLE = config['sample']
cores = config.get("max_cores", 16)  # Default to 16 if not set


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
        "../envs/methylcibersort.yaml"
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

