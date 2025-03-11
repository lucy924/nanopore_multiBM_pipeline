SAMPLE = config['sample']

rule combine_bedmethyls:
    input:
        bedmeth1 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.1.bedmethyl.gz",
        bedmeth2 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.2.bedmethyl.gz",
        bedmeth3 = f"results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_mods.ungrouped.bedmethyl.gz"
    output:
        temp(f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed")
    params:
        threads = config['threads']
    shell:
        "bgzip -dc {input.bedmeth1} {input.bedmeth2} {input.bedmeth3} | "
        "sort --parallel={params.threads} -k1,1 -k2,2n > {output}"

rule convert_bedmethyl_to_DSS:
    input: 
        f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.bedmethyl.bed"
    output: 
        f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    script:
        "scripts/convert_DSS.py"


rule transform_DSS_for_ImIn:
    input:
        dss_file = f"results/{SAMPLE}/mod_calling/{SAMPLE}.wf_mods.all.dss_format.tsv"
    params:
        sample_name = SAMPLE

# """
# awk -v OFS='\t' 'BEGIN{print "chr","pos","N","X"}{print $1,$2,($12+$13),$13}' {input} > {output}
# """
