#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# sv_annotation.py
################################################################################


import pandas as pd
from cyvcf2 import VCF
from snakemake.script import snakemake
from shared_functions import variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict, get_snp_by_genomic_location, variant_dict_columns_to_add

# -------------------------------------------------------------------------- #
# Data from snakefile
panel_metadata_fp = snakemake.input['panel_metadata']
vcf_sv_fp = snakemake.input['vcf_sv_gz']
sv_output = snakemake.output['sv_csv']

log = open(snakemake.log[0], 'w')
# ------------------------------------------------ #

variants_metadata_df_svs = variant_prep(panel_metadata_fp, variant_type = 'sv')

vcf_sv = VCF(vcf_sv_fp)

info_to_add_to_metadata = dict()
num_variants_found_total = 0
rows_found = list()
rows_not_found = list()

for j, row in enumerate(variants_metadata_df_svs.iterrows()):
    target_ID = row[1]['ID']
    if target_ID in rows_found:
        continue
    else:
        entry_found = False
    
    print(f"target {target_ID} of snp metadata")
    
    if target_ID not in info_to_add_to_metadata.keys():
        info_to_add_to_metadata[target_ID] = dict()
        # Initialise columns so they are added regardless of the result
        for colname in variant_dict_columns_to_add:
            info_to_add_to_metadata[target_ID][colname] = ''
        
    # Get the location info as a string
    loc, chrom, start, end = get_location_string(row_data = row[1])
    loc_list = [chrom, start, end]
    
    in_vcf = False  # for checking if the location is in the vcf but entry isn't found
    for i, variant in enumerate(vcf_sv(loc)):
        print('i: ', i)
        in_vcf = True
        info_to_add_to_metadata[target_ID]['ClinVar'] = 'N/A'

        if len(variant.genotypes) > 1:
            raise ValueError("Multiple sample processing not currently supported.")
        
        # Get info into dict so I can work with it
        info_dict = get_annotation_info_dict(info_field = variant.INFO)
        annotation_dict = get_annotation_dict(info_dict_ann = info_dict['ANN'])

        # 1. Match genomic location
        if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
            info_to_add_to_metadata[target_ID], entry_found = get_snp_by_genomic_location(
                variants_metadata_df_svs, variant, info_dict, annotation_dict, 
                info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list
            )
        if entry_found:
            rows_found.append(target_ID)
            break
        
        # 3. Match end genomic position
        if (chrom == variant.CHROM) and (end == variant.end):
            info_to_add_to_metadata[target_ID], entry_found = get_snp_by_genomic_location(
                variants_metadata_df_svs, variant, info_dict, annotation_dict, 
                info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list, end_only = True
            )
        if entry_found:
            rows_found.append(target_ID)
            break
            
        print('================')

    if entry_found:
        # move out into next row of metadata
        num_variants_found_total += 1
        rows_found.append(target_ID)
        continue
    elif in_vcf:
        print("area is in vcf but entry not found:")
        print(row[0])
        print(row[1])
        raise ValueError()
    else:
        rows_not_found.append(target_ID)


log.write(f"number of variants found total: {num_variants_found_total}\n")
log.write(f"number of variants not found: {len(rows_not_found)}\n")
log.write(f'number of variants I was looking for: {len(variants_metadata_df_svs)}\n')
print(f"number of variants found total: {num_variants_found_total}")
print(f"number of variants not found: {len(rows_not_found)}")
print(f'number of variants I was looking for: {len(variants_metadata_df_svs)}')

# Convert dictionary to DataFrame
info_df = pd.DataFrame(info_to_add_to_metadata).T
info_df.index.name = 'ID'
info_df.reset_index(inplace=True)

# Merge the DataFrames
merged_df = pd.merge(variants_metadata_df_svs, info_df, on='ID', how='left')

# Output snp df to csv
merged_df.to_csv(sv_output, index = False)

log.close()
