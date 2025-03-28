#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# shared_functions.py
################################################################################

import pandas as pd

HG_LENGTH = 3100000000
# HG_length = 3000000000

CHROMOSOMES = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
]

CHROM_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569
}

VARIANT_TYPE = "Variant Type (snp, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, demographic, clinicopathology)"

preclin_stage_panel_result_header = ["ID", "Marker name", "Scoring Type", "Biomarker Type", "Result Options", "Result"]

variant_dict_columns_to_add = ['ClinVar', 'Significance (ClinVar)', 'Consequence (Clinvar)', 'Reference Allele', 'Variant Allele', 'Genotype', 'HGVS.c', 'HGVS.p', 'SV Length', 'SV Type']

def filter_to_variant_type(panel_metadata_df, variant_type):
    """
    Filter rows based on Variant type entry.
    Variant type must be one of: snp, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, microsatellite, demographic, clinicopathology
    """
    # Sanity check:
    variant_types = ["snp", "sv", "mod", "area_mutations", "expression", "exp_ratio", "immune_ratio", "immune_inf", "microsatellite", "demographic", "clinicopathology"]
    if variant_type not in variant_types:
        raise ValueError("variant_type must be one of: snp, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, microsatellite, demographic, clinicopathology")
    
    panel_metadata_df = panel_metadata_df[
        (panel_metadata_df[VARIANT_TYPE] == variant_type)
    ]
    
    #TODO - for my sanity during development
    # if variant_type == 'snp':
    #     panel_metadata_df = panel_metadata_df[
    #         (panel_metadata_df["length"] == '1')
    #     ]  # will need to deal with these caveats down the track
        
    return panel_metadata_df


def reduce_metadata_df(df):
    """Drop columns of excess info from metadata for snvs and svs"""
    df = df.drop(
        columns=[
            # "Is this record the whole gene?",
            "Is variant in coding region?",
            "Illumina EPIC ID",
            # "DNA methylation",
            "Notes",
            "WARNINGS",
            "Ref",
            "length",
        ]
    )
    return df


def variant_prep(path2_variants_metadata_csv, variant_type):
    """
    Load in metadata files. 
    """
    import pandas as pd
    
    with open(path2_variants_metadata_csv, "r") as fpc:
        variants_metadata_df = pd.read_csv(fpc, dtype={"ID": str})

    # clean up csv of empty columns
    columns_to_drop = variants_metadata_df.filter(like="Unnamed").columns
    variants_metadata_df.drop(columns=columns_to_drop, inplace=True)
    
    variants_metadata_df_v = filter_to_variant_type(variants_metadata_df, variant_type = variant_type)
    variants_metadata_df_v = reduce_metadata_df(variants_metadata_df_v)
    
    return variants_metadata_df_v


def get_location_string(row_data):
    chrom = row_data['chrom']
    start = int(row_data['start pos'].replace(",", ""))
    end = int(row_data['end pos'].replace(",", ""))
    loc = (f'{chrom}:{start - 2}-{end + 2}')
    return(loc, chrom, start, end)


def get_annotation_dict(info_dict_ann):
    # from vcf header info lines
    annotation_header = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
    annotation_dict = dict(zip(annotation_header.split(' | '), info_dict_ann.split('|')))
    return(annotation_dict)
    
    
def get_annotation_info_dict(info_field):
    info_dict = {}
    for entry in info_field:
        info_dict[entry[0]] = entry[1]
        
    return info_dict


def add_result(add_data_to_this_dict, variant, annotation_dict, log, clinvar=False):
    from warnings import warn
    
    if not variant.FILTERS[0] == 'PASS':
        log.write(f"  WARNING: This variant did not pass whatever filters there are in the VCF.\nVariant position: {variant.CHROM}:{variant.start}-{variant.end}")
        warn(f"This variant did not pass whatever filters there are in the VCF.\nVariant position: {variant.CHROM}:{variant.start}-{variant.end}")
    
    ref_allele = variant.REF
    
    # Extract and print the genotype information
    if len(variant.genotypes) > 1:
        log.write("Multiple sample processing not currently supported.\n")
        raise ValueError("Multiple sample processing not currently supported.")
    genotype_access_list = [variant.REF]
    genotype_access_list.extend(variant.ALT)
    for genotype in variant.genotypes:
        allele1, allele2, phased = genotype
        allele_base_1 = genotype_access_list[allele1]
        allele_base_2 = genotype_access_list[allele2]
        genotype_text = (f"{allele_base_1}|{allele_base_2}")
    
    try:
        annotation_consequence = annotation_dict['Annotation']  # missense_variant
    except:
        annotation_consequence = 'N/Av'
    try:
        significance = variant.INFO.get('CLNSIG')
    except:
        significance = 'N/Av'
    try:
        sv_length = variant.INFO.get('SVLEN') 
    except:
        sv_length = 'N/A'
    try:
        sv_type = variant.INFO.get('SVTYPE') 
    except:
        sv_type = 'N/A'

    annotation_HGVSc = annotation_dict['HGVS.c']
    annotation_HGVSp = annotation_dict['HGVS.p']
        
    add_data_to_this_dict['Significance (ClinVar)'] = significance
    add_data_to_this_dict['Consequence (Clinvar)'] = annotation_consequence
    add_data_to_this_dict['Reference Allele'] = ref_allele
    add_data_to_this_dict['Genotype'] = genotype_text
    add_data_to_this_dict['HGVS.c'] = annotation_HGVSc
    add_data_to_this_dict['HGVS.p'] = annotation_HGVSp
    add_data_to_this_dict['SV Length'] = sv_length
    add_data_to_this_dict['SV Type'] = sv_type

    return add_data_to_this_dict


def get_snp_by_genomic_location(variants_metadata_df_sub, variant, info_dict, annotation_dict, row_data, entry_found, log, location = list(), end_only=False, start_only=False, clinvar = False):
    """
    getting snp by full or just 'start'/'end' genomic location.  
    I think if start only, likely an indel
    location = [chrom, start, end]
    """
    [chrom, start, end] = location
    
    if end_only:
        filterdata_loc = variants_metadata_df_sub.loc[
            (variants_metadata_df_sub['chrom'] == chrom) &
            (variants_metadata_df_sub['end pos'] == "{:,}".format(end))  # formatted to string
        ]    
    elif start_only:
        filterdata_loc = variants_metadata_df_sub.loc[
            (variants_metadata_df_sub['chrom'] == chrom) &
            (variants_metadata_df_sub['start pos'] == "{:,}".format(start))  # formatted to string
        ]
    else:
        filterdata_loc = variants_metadata_df_sub.loc[
            (variants_metadata_df_sub['chrom'] == chrom) &
            (variants_metadata_df_sub['start pos'] == "{:,}".format(start)) &
            (variants_metadata_df_sub['end pos'] == "{:,}".format(end))
        ]
    
    if len(filterdata_loc) == 1:
        if end_only:
            log.write('  entry found by end only genomic location\n')
        elif start_only:
            log.write('  entry found by start only genomic location - maybe an indel\n')
        else:
            log.write('  entry found by genomic location\n')
        
        # Actually do the thing
        row_data = add_result(
            add_data_to_this_dict = row_data, variant = variant, 
            annotation_dict = annotation_dict, log = log, clinvar = clinvar)
        entry_found = True
    
    else:
        log.write(
            f"Location matches but won't give me a single entry for some reason.\nVariant position: {chrom}:{start}-{end}\n")
        raise ValueError(
            f"Location matches but won't give me a single entry for some reason.\nVariant position: {chrom}:{start}-{end}")
    
    return row_data, entry_found


def add_functional_flanking_regions(panel_csv):
    promoter_length = 2000
    end_length = 1000

    for index, entry in panel_csv.iterrows():
        # only look at whole genes to add promotor and end regions
        if entry['Is this record the whole gene?'] == 'Yes':
            # ensure positions are integers
            if isinstance(entry['start pos'], str) and ',' in entry['start pos']:
                start_coord = int(entry['start pos'].replace(',', ''))
            else:
                start_coord = int(entry['start pos'])
            if isinstance(entry['end pos'], str) and ',' in entry['end pos']:
                end_coord = int(entry['end pos'].replace(',', ''))
            else:
                end_coord = int(entry['end pos'])
            
            new_start_coord = None
            new_end_coord = None
            strand = entry['strand']
            
            if strand == '+':
                gene_start_coord = start_coord
                gene_end_coord = end_coord
                new_start_coord = gene_start_coord - promoter_length
                new_end_coord = gene_end_coord + end_length
            elif strand == '-':
                gene_start_coord = end_coord
                gene_end_coord = start_coord
                new_start_coord = gene_end_coord - end_length
                new_end_coord = gene_start_coord + promoter_length
            
            # update coords
            panel_csv.at[index, 'start pos'] = new_start_coord
            panel_csv.at[index, 'end pos'] = new_end_coord
        
        # replace empty strand with *
        if entry['strand'] not in ["+", "-"]:
            panel_csv.at[index, 'strand'] = "*"
            
    return panel_csv
