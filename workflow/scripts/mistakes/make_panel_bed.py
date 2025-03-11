#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# make_panel_bed.py
################################################################################

from warnings import WarningMessage
import numpy as np
import pandas as pd
from snakemake.script import snakemake


def add_functional_flanking_regions(panel_csv):
    promoter_length = 2000
    end_length = 1000
    default_gene_buff = 2000

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
            else:
                # replace empty strand with *
                panel_csv.at[index, 'strand'] = "*"
                new_start_coord = gene_start_coord - default_gene_buff
                new_end_coord = gene_end_coord + default_gene_buff
            
            # update coords
            panel_csv.at[index, 'start pos'] = new_start_coord
            panel_csv.at[index, 'end pos'] = new_end_coord
            
    return panel_csv

def add_buff_regions_to_snvs(panel_csv):
    default_buff = 50
    for index, entry in panel_csv.iterrows():
        # only look at whole genes to add promotor and end regions
        if entry['Is this record the whole gene?'] != 'Yes':
            # ensure positions are integers
            if isinstance(entry['start pos'], str) and ',' in entry['start pos']:
                start_coord = int(entry['start pos'].replace(',', ''))
            else:
                start_coord = int(entry['start pos'])
            if isinstance(entry['end pos'], str) and ',' in entry['end pos']:
                end_coord = int(entry['end pos'].replace(',', ''))
            else:
                end_coord = int(entry['end pos'])
            
            # Add default buff for all other entries
            new_start_coord = start_coord - default_buff
            new_end_coord = end_coord + default_buff
            
        # update coords
            panel_csv.at[index, 'start pos'] = new_start_coord
            panel_csv.at[index, 'end pos'] = new_end_coord
            
    return panel_csv

def restructure_to_bed(df):
    """ add empty score column and reorder columns"""
    df['score'] = np.nan
    df = df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
    return df


def get_immune_infiltrate_ref_locs(immune_reference_dataset, epic_locs_hg38):    
    
    # Load immune reference dataset
    with open(immune_reference_dataset, 'r') as fp:
        imm_ref = pd.read_csv(fp)
    probe_list = imm_ref['CpGs']
    
    # Load epic genome locations
    with open(epic_locs_hg38, 'r') as fp:
        EPIC_probes_loc_hg38 = pd.read_csv(fp, index_col=0)
        
    probes_genomic_locs = EPIC_probes_loc_hg38[EPIC_probes_loc_hg38['probe'].isin(probe_list)]
    
    if len(probes_genomic_locs) < len(probe_list):
        WarningMessage(
            f'not all genomic locations were found for reference probes ({len(probe_list) - len(probes_genomic_locs)} are missing.)',
            category = UserWarning,
            filename = 'make_panel_bed.py',
            lineno = 71)
    
    # rename existing columns
    immune_probes_bed_df = probes_genomic_locs.rename(columns = {
        'probe': 'name',
        'seqnames': '#chrom',
        'start': 'chromStart',
        'end': 'chromEnd'
    })
    immune_probes_bed_df = restructure_to_bed(immune_probes_bed_df)
    
    return immune_probes_bed_df


def add_immune_infiltrate_locations(input_bed, immune_reference_dataset, epic_locs_hg38):
    
    immune_probes_bed_df = get_immune_infiltrate_ref_locs(immune_reference_dataset, epic_locs_hg38)
    
    # Merge variants and probes bed files
    merged_df = pd.concat([input_bed, immune_probes_bed_df])

    return merged_df


# ------------------------------------------------ #
# Data from snakefile
panel_csv_fp = snakemake.input['panel_csv']
immune_reference_dataset = snakemake.input['immune_reference_dataset']
epic_locs_hg38 = snakemake.input['epic_locs_hg38']
# panel_genebuffed_csv_fp = snakemake.output['panel_genebuffed_csv']
panel_bed_fp = snakemake.output['panel_bed']
all_targets_fp = snakemake.output['all_targets']

# ------------------------------------------------ #
# Load panel
with open(panel_csv_fp, 'r') as fp:
    panel_csv = pd.read_csv(fp, dtype={'ID': str}, thousands = ',')

# ------------------------------------------------ #
# Add functional flanking regions
panel_csv = add_functional_flanking_regions(panel_csv)
panel_csv.to_csv(snakemake.output['panel_csv_b4_targets'])
panel_csv_for_targets_file = add_buff_regions_to_snvs(panel_csv.copy())
panel_csv_for_targets_file.to_csv(snakemake.output['panel_csv_targets'])

# ------------------------------------------------ #
# reformat panels
panel_bed = panel_csv[['ID', 'chrom', 'start pos', 'end pos', 'strand']]
panel_bed = panel_bed.rename(columns = {
    'ID': 'name', 
    'chrom': '#chrom', 
    'start pos': 'chromStart', 
    'end pos':'chromEnd'})
panel_bed = restructure_to_bed(panel_bed)

panel_bed_for_targets = panel_csv_for_targets_file[['ID', 'chrom', 'start pos', 'end pos', 'strand']]
panel_bed_for_targets = panel_bed_for_targets.rename(columns = {
    'ID': 'name', 
    'chrom': '#chrom', 
    'start pos': 'chromStart', 
    'end pos':'chromEnd'})
panel_bed_for_targets = restructure_to_bed(panel_bed_for_targets)

# ------------------------------------------------ #
# Save panel for post-seq analysis
# panel_csv.to_csv(panel_genebuffed_csv_fp, index = False)  # only genes are buffed
panel_bed.to_csv(panel_bed_fp, sep = '\t', index = False)  # only genes are buffed

# ------------------------------------------------ #
# Add immune infiltrate probes to bed for minknow file
all_targets = add_immune_infiltrate_locations(panel_bed_for_targets, immune_reference_dataset, epic_locs_hg38)

# ------------------------------------------------ #
# Save targets for adding sequence buffer regions
all_targets.to_csv(all_targets_fp, sep = '\t', index = False)
