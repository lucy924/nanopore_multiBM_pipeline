
SAMPLE = config['sample']
cores = config.get("max_cores", 16)  # Default to 16 if not set


rule bgzip_clinvar_vcf:
    input:
        vcf_clinvar = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf"
    output:
        vcf_clinvar_gz = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz"
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.bgzip_clinvar_vcf.tsv"
    shell:
        "bgzip -k {input}"


rule tabix_clinvar_vcf_gz:
    input:
        vcf_clinvar_gz = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz"
    output:
        vcf_clinvar_gz_tbi = f"results/{PROJECT}/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz.tbi"
    benchmark:
        f"logs/{PROJECT}/{SAMPLE}/{SAMPLE}.benchmark.tabix_clinvar_vcf_gz.tsv"
    shell:
        "tabix {input}"

