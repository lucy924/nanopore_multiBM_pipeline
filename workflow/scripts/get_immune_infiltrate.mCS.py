#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# get_immune_infiltrate.py
################################################################################


import os
import json
import numpy as np
import pandas as pd
from snakemake.script import snakemake  # type: ignore
from shared_functions import VARIANT_TYPE, preclin_stage_panel_result_header, variant_prep, BIOMARKER_NAME, SCORING_TYPE, RESULT_OPTIONS

lymphocytes = ["CD19", "CD4_Eff",
               "CD56", "CD8", "Treg", ]
monocytes = ["CD14"]
neutrophils = ["Neu"]
eosinophils = ["Eos"]


def get_LMR(deconv_df, log):
    """Lymphocyte/Monocyte ratio"""
    lymphocyte_val = sum(deconv_df[lymphocytes].iloc[0])
    monocyte_val = sum(deconv_df[monocytes].iloc[0])

    # Override for testing
    # monocyte_val = 0.015

    if monocyte_val < 0.0001:
        log.write("Monocyte value too low to get LMR.\n")
        LMR_ratio = np.nan
        
    elif lymphocyte_val < 0.0001:
        log.write("Lymphocyte value too low to get LMR.\n")
        LMR_ratio = np.nan
        
    else:
        LMR_ratio = lymphocyte_val / monocyte_val
        log.write(f"Lymphocyte to Monocyte ratio (LMR) is {LMR_ratio:.2f}\n")
        
    return LMR_ratio


def get_NLR(deconv_df, log):
    """Neutrophil/Lymphocyte ratio"""
    lymphocyte_val = sum(deconv_df[lymphocytes].iloc[0])
    neutrophil_val = sum(deconv_df[neutrophils].iloc[0])
    
    if neutrophil_val < 0.0001:
        log.write("Neutrophil value too low to get NLR.\n")
        NLR_ratio = np.nan
    
    elif lymphocyte_val < 0.0001:
        log.write("Lymphocyte value too low to get NLR.\n")
        NLR_ratio = np.nan
        
    else:
        NLR_ratio = neutrophil_val / lymphocyte_val
        log.write(f"Neutrophil to Lymphocyte ratio (NLR) is {NLR_ratio:.2f}\n")
    
    return NLR_ratio


def get_dendritic_cells():
    return


def get_Th1_cells():
    return


def get_M2_macrophages():
    return


def get_platelet_count():
    return

# ------------------------------------------------ #
# get snakemake variables
deconv_fp = snakemake.input["mCS_results"]
immune_results_fp = snakemake.output["immune_results"]

# path2_panel_metadata_csv = snakemake.input["panel_metadata"]
# results_fp = snakemake.output['panel_mod_results']

log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #

# Load docker CS output
with open(deconv_fp, "r") as fd:
    deconv_df = pd.read_csv(fd, sep='\t')

sample_row_id = 0

Monocytes = deconv_df["CD14"][sample_row_id]
Bcells = deconv_df["CD19"][sample_row_id]
CD4_Tcells = deconv_df["CD4_Eff"][sample_row_id]
NK_cells = deconv_df["CD56"][sample_row_id]
CD8_Tcells = deconv_df["CD8"][sample_row_id]
Tregs = deconv_df["Treg"][sample_row_id]
Endothelial = deconv_df["Endothelial"][sample_row_id]
Eosinophils = deconv_df["Eos"][sample_row_id]
Fibroblasts = deconv_df["Fibroblast"][sample_row_id]
Neutrophils = deconv_df["Neu"][sample_row_id]
Cancer = deconv_df["Cancer"][sample_row_id]

# Get ratios of interest
LMR_ratio = get_LMR(deconv_df, log)
NLR_ratio = get_NLR(deconv_df, log)

log.write(f'Proportion of infiltrating Monocytes is {Monocytes}\n')
log.write(f'Proportion of infiltrating B-cells is {Bcells}\n')
log.write(f'Proportion of infiltrating CD4+ T-cells is {CD4_Tcells}\n')
log.write(f'Proportion of infiltrating CD8+ T-cells is {NK_cells}\n')
log.write(f'Proportion of infiltrating NK cells is {CD8_Tcells}\n')
log.write(f'Proportion of infiltrating Neutrophils is {Neutrophils}\n')
log.write(f'Proportion of infiltrating Tregs is {Tregs}\n')
log.write(f'Proportion of infiltrating Endothelial cells is {Endothelial}\n')
log.write(f'Proportion of infiltrating Eosinophils is {Eosinophils}\n')
log.write(f'Proportion of infiltrating Fibroblasts is {Fibroblasts}\n')
log.write(f'Proportion of Cancer cells is {Cancer}\n')
log.write(f'')
log.write(f'Lymphocyte to Monocyte ratio is: {LMR_ratio}\n')
log.write(f'Neutrophil to Lymphocyte ratio is {NLR_ratio}\n')

# Load panel info
panel_data_ratio = variant_prep(snakemake.input["panel_metadata"], 'immune_ratio')
panel_data_infiltrate = variant_prep(snakemake.input["panel_metadata"], 'immune_inf')

# Make preclin panel output
preclin_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)

for i, row in panel_data_ratio.iterrows():
    if row[BIOMARKER_NAME] == "LMR":
        result = LMR_ratio
    elif row[BIOMARKER_NAME] == 'NLR':
        result = NLR_ratio
    else:
        result = np.nan
        
    preclin_panel_df.loc[i] = [row['ID'], row[BIOMARKER_NAME], row[SCORING_TYPE], row[VARIANT_TYPE], row[RESULT_OPTIONS], result]  
    
for i, row in panel_data_infiltrate.iterrows():
    
    if row[BIOMARKER_NAME] == 'Monocyte_inf':
        result = Monocytes
    elif row[BIOMARKER_NAME] == 'Bcell_inf':
        result = Bcells
    elif row[BIOMARKER_NAME] == 'CD4_inf':
        result = CD4_Tcells
    elif row[BIOMARKER_NAME] == 'NK_inf':
        result = NK_cells
    elif row[BIOMARKER_NAME] == 'CD8_inf':
        result = CD8_Tcells
    elif row[BIOMARKER_NAME] == 'Treg_inf':
        result = Tregs
    elif row[BIOMARKER_NAME] == 'Neutrophil_inf':
        result = Neutrophils
    elif row[BIOMARKER_NAME] == 'Endothelial_inf':
        result = Endothelial
    elif row[BIOMARKER_NAME] == 'Eosinophil_inf':
        result = Eosinophils
    elif row[BIOMARKER_NAME] == 'Fibroblast_inf':
        result = Fibroblasts
    elif row[BIOMARKER_NAME] == 'Cancer_inf':
        result = Cancer
    else:
        result = np.nan
        
    preclin_panel_df.loc[i] = [row['ID'], row[BIOMARKER_NAME], row[SCORING_TYPE], row[VARIANT_TYPE], row[RESULT_OPTIONS], result]  

# save preclin panel
preclin_panel_df.to_csv(immune_results_fp, index=False)

log.close()

# Notes:

# Platelet isn't available as they don't have gDNA. Future investigations may reveal the ability to use mDNA signatures for this
# https://www.sciencedirect.com/science/article/pii/S0167527323016960
# Th1 isn't available with this reference
# Dendritic cells aren't available with this reference
# M2 macrophages specifically aren't available with this reference
