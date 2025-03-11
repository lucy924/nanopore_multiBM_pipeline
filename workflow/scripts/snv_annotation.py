#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# sv_annotation.py
################################################################################

import os
if os.getenv("SNAKEMAKE_DEBUG"):
    class FakeSnakemake:
        SAMPLE = "test3"
        PREFIX = "/external/analyses/lucy/nanopore_multiBM_pipeline"
        input = {
            'panel_metadata': f"{PREFIX}/config/testing_panel_metadata.csv",
            'vcf_clinvar_gz': f"{PREFIX}/results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz",
            'vcf_clinvar_gz_tbi': f"{PREFIX}/results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp_clinvar.vcf.gz.tbi",
            'vcf_all': f"{PREFIX}/results/{SAMPLE}/wf-humvar/{SAMPLE}.wf_snp.vcf.gz",
            }
        output = {
            'snv_csv': f"{PREFIX}/results_debug/{SAMPLE}/snv_annotation/{SAMPLE}.raw_snv_results.csv"
            }
        log = [f"{PREFIX}/results_debug/{SAMPLE}.snv_annotation.log"]

    snakemake = FakeSnakemake()
    
from warnings import warn
import pandas as pd
from cyvcf2 import VCF
# from snakemake.script import snakemake
from shared_functions import variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict, add_result, get_snp_by_genomic_location, variant_dict_columns_to_add, VARIANT_TYPE,preclin_stage_panel_result_header


def get_snp_by_RS_number(variants_metadata_df_sub, variant, info_dict, annotation_dict, row_data, entry_found, log, clinvar = False):
    """
    VCF header:
    ##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
    RS numbers have the same position entry for start and end, and it is the end number. 
    See help pic above.
    """
    
    dbSNP_ID = info_dict['RS']
    
    filterdata_dbSNPid = variants_metadata_df_sub.loc[variants_metadata_df_sub['SNP ID'] == f'rs{dbSNP_ID}']
    
    if len(filterdata_dbSNPid) == 1:
        log.write('  entry found by dbSNPid\n')
        row_data = add_result(
            add_data_to_this_dict = row_data, variant = variant, annotation_dict = annotation_dict, log = log, clinvar = clinvar)
        entry_found = True
    else:
        log.write('snpid not found when it should be (line 86)')
        raise ValueError('snpid not found when it should be (line 86)')
        
    return row_data, entry_found


# -------------------------------------------------------------------------- #
# Data from snakefile
panel_metadata_fp = snakemake.input['panel_metadata']
vcf_clinvar_fp = snakemake.input['vcf_clinvar_gz']
vcf_all_fp = snakemake.input['vcf_all']
snp_output = snakemake.output['snv_csv']
snp_preclin_output = snakemake.output['snv_panel_csv']

log = open(snakemake.log[0], 'w')
# ------------------------------------------------ #
variants_metadata_df_snps = variant_prep(panel_metadata_fp, 'snp')

vcf_clinvar = VCF(vcf_clinvar_fp)
vcf_all = VCF(vcf_all_fp)
vcfs = [vcf_clinvar, vcf_all]

info_to_add_to_metadata = dict()
num_variants_found_total = 0
num_variants_found_clinvar = 0
rows_found = list()
rows_not_found = list()

breaking = False
for i, vcf in enumerate(vcfs):
    if i == 0:
        clinvar = True
        log.write('################## CLINVAR ######################\n')
    else:
        clinvar = False
        log.write('################## NOT CLINVAR ######################\n')
        
    for j, row in enumerate(variants_metadata_df_snps.iterrows()):
        target_ID = row[1]['ID']
        if target_ID in rows_found:
            continue
        else:
            entry_found = False
        
        log.write(f"looking for panel target {target_ID} of snp metadata\n")
        
        if target_ID not in info_to_add_to_metadata.keys():
            info_to_add_to_metadata[target_ID] = dict()
            # Initialise columns so they are added regardless of the result
            for colname in variant_dict_columns_to_add:
                info_to_add_to_metadata[target_ID][colname] = ''

        # Get the location info as a string
        loc, chrom, start, end = get_location_string(row_data = row[1])
        loc_list = [chrom, start, end]
        
        for variant in vcf(loc):
            add_variant_to_data = ()
            log.write('================\n')
            if clinvar: 
                info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes'
                num_variants_found_clinvar += 1    
            else:
                # Can add here because row would have been skipped if it was found in clinvar
                info_to_add_to_metadata[target_ID]['ClinVar'] = 'No'

            if len(variant.genotypes) > 1:
                log.write("Multiple sample processing not currently supported.\n")
                raise NotImplementedError("Multiple sample processing not currently supported.")
                            
            # Get info into dict so I can work with it
            info_dict = get_annotation_info_dict(info_field = variant.INFO)
            annotation_dict = get_annotation_dict(info_dict_ann = info_dict['ANN'])

            # 1. Matching RS number
            if 'RS' in info_dict.keys():
                info_to_add_to_metadata[target_ID], entry_found = get_snp_by_RS_number(
                    variants_metadata_df_snps, variant, info_dict, annotation_dict, 
                    info_to_add_to_metadata[target_ID], entry_found, log, clinvar = clinvar
                )
            if entry_found:
                rows_found.append(target_ID)
                break
                
            # 2. Match genomic location
            if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
                info_to_add_to_metadata[target_ID], entry_found = get_snp_by_genomic_location(
                    variants_metadata_df_snps, variant, info_dict, annotation_dict, 
                    info_to_add_to_metadata[target_ID], entry_found, log = log, location = loc_list,
                    clinvar = clinvar
                )
            if entry_found:
                rows_found.append(target_ID)
                break
            
            # 3. Match end genomic position only
            if (chrom == variant.CHROM) and (end == variant.end):
                info_to_add_to_metadata[target_ID], entry_found = get_snp_by_genomic_location(
                    variants_metadata_df_snps, variant, info_dict, annotation_dict, info_to_add_to_metadata[target_ID], entry_found, log, location = loc_list, end_only = True, 
                    clinvar = clinvar
                )
            if entry_found:
                rows_found.append(target_ID)
                break
            
            # 4. Match start genomic position only - likely indel
            if (chrom == variant.CHROM) and (start == variant.start):
                info_to_add_to_metadata[target_ID], entry_found = get_snp_by_genomic_location(
                    variants_metadata_df_snps, variant, info_dict, annotation_dict, info_to_add_to_metadata[target_ID], entry_found, log, location = loc_list, end_only = False, start_only = True,
                    clinvar = clinvar
                )
            if entry_found:
                rows_found.append(target_ID)
                break
            else:
                log.write("in vcf but not found yet:\n")
                log.write(str(row[0]) + '\n')
                log.write(str(row[1]) + '\n')
                print("in vcf but not found yet:\n")
                print(row[0])
                print(row[1])
                raise ValueError()
        log.write('\n')

        if entry_found:
            # move out into next row of metadata
            num_variants_found_total += 1
            rows_found.append(target_ID)
            continue
        else:
            if not clinvar:  
                # don't need to do it with clinvar as it has not had a chance to search the normal vcf
                rows_not_found.append(target_ID)

log.write(f"\nnumber of variants found total: {num_variants_found_total}\n")
log.write(f"  number of variants found in clinvar: {num_variants_found_clinvar}\n")
log.write(f"number of variants not found: {len(rows_not_found)}\n")
log.write(f'number of variants I was looking for: {len(variants_metadata_df_snps)}\n')
print(f"number of variants found total: {num_variants_found_total}")
print(f"  number of variants found in clinvar: {num_variants_found_clinvar}")
print(f"number of variants not found: {len(rows_not_found)}")
print(f'number of variants I was looking for: {len(variants_metadata_df_snps)}')

# Convert dictionary to DataFrame
info_df = pd.DataFrame(info_to_add_to_metadata).T
info_df.index.name = 'ID'
info_df.reset_index(inplace=True)

# Merge the DataFrames
merged_df = pd.merge(variants_metadata_df_snps, info_df, on='ID', how='left')

# Output snp df to csv
merged_df.to_csv(snp_output, index = False)

# Get preclin_panel output
preclin_stage_panel_result_header = ["ID", "Scoring Type", "Biomarker Type", "Result Options", "Result"]
preclin_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)
# Set the dtypes for the columns
preclin_panel_df = preclin_panel_df.astype(str)

# Only get panel id with stuff in genotype result
only_genotypes = merged_df[pd.notna(merged_df['Genotype'])]
for i, row in only_genotypes.iterrows():
    preclin_panel_df.loc[i] = [row['ID'], row['Scoring Type'], row[VARIANT_TYPE], row['Variant'], row['Genotype']]

preclin_panel_df.to_csv(snp_preclin_output, index = False)

log.close()
