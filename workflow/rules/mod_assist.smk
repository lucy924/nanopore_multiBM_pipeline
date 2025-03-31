
SAMPLE = config['sample']
cores = config.get("max_cores", 16)  # Default to 16 if not set


rule combine_bedmethyls:
    input:
        bedmeth1 = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.1.bedmethyl.gz",
        bedmeth2 = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.2.bedmethyl.gz",
        bedmeth3 = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.ungrouped.bedmethyl.gz"
    output:
        f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed"
    params:
        threads = cores
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.combine_bedmethyls.tsv"
    shell:
        (
        "bgzip -dc {input.bedmeth1} {input.bedmeth2} {input.bedmeth3} | "
        "sort --parallel={params.threads} -k1,1 -k2,2n > {output}"
        )


rule convert_bedmethyl_to_DSS:
    input: 
        f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed"
    output: 
        f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.convert_bedmethyl_to_DSS.tsv"
    script:
        "../scripts/convert_DSS.py"


rule prep_for_getting_betas:
    input:
        IlluminaEPIC_genomic_locations_hg38 = "resources/IlluminaEPIC_genomic_locations_hg38.csv",
        dss_file = f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    output:
        pre_beta = f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.pre_beta.csv"
    log:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.prep_for_betas.log"
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.prep_for_getting_betas.tsv"
    script:
        "../scripts/process_result_before_betas.py"

rule add_betas:
    input:
        pre_beta = f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.pre_beta.csv"
    output:
        # post_beta = temp(f"results/{SAMPLE}/mod_calling/{SAMPLE}.post_beta.csv")
        post_beta = f"results/{PROJECT}/{SAMPLE}/mod_calling/{SAMPLE}.post_beta.csv"
    log:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.add_betas.log"
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.add_betas.tsv"
    conda: 
        "../envs/methylcibersort.yaml"
    script:
        "../scripts/add_betas.R"

