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
from snakemake.script import snakemake
from shared_functions import CHROMOSOMES

# lymphocytes = ["B-cells_EPIC", "CD4T-cells_EPIC",
#                "NK-cells_EPIC", "CD8T-cells_EPIC"]
# monocytes = ["Monocytes_EPIC"]
# neutrophils = ["Neutrophils_EPIC"]

lymphocytes = ["CD19", "CD4_Eff",
               "CD56", "CD8", "Treg", ]
monocytes = ["CD14"]
neutrophils = ["Neu"]
eosinophils = ["Eos"]


def get_LMR(deconv_df, log):
    """Lymphocyte/Monocyte ratio"""
    lymphocyte_val = deconv_df[deconv_df["Immune cells"].isin(
        lymphocytes)]["sample % cell"].sum()
    
    monocyte_val = deconv_df[deconv_df["Immune cells"].isin(
        monocytes)]["sample % cell"].sum()

    # Override for testing
    # monocyte_val = 0.015

    if monocyte_val < 0.001:
        log.write("Monocyte value too low to get LMR.\n")
        LMR_ratio = np.nan
        
    elif lymphocyte_val < 0.001:
        log.write("Lymphocyte value too low to get LMR.\n")
        LMR_ratio = np.nan
        
    else:
        LMR_ratio = lymphocyte_val / monocyte_val
        log.write(f"Lymphocyte to Monocyte ratio (LMR) is {LMR_ratio:.2f}\n")
        
    return LMR_ratio


def get_NLR(deconv_df, log):
    """Neutrophil/Lymphocyte ratio"""
    lymphocyte_val = deconv_df[deconv_df["Immune cells"].isin(
        lymphocytes)]["sample % cell"].sum()
    
    neutrophil_val = deconv_df[deconv_df["Immune cells"].isin(
        neutrophils)]["sample % cell"].sum()
    
    if neutrophil_val < 0.001:
        log.write("Neutrophil value too low to get NLR.\n")
        NLR_ratio = np.nan
    
    elif lymphocyte_val < 0.001:
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
immune_reference_dataset = snakemake.input['immune_reference_dataset']
meth_results_fp = snakemake.input["epic_probe_results"]
deconv_fp = snakemake.input["methatlas_deconvolution"]
results_json_imin = snakemake.output["results_json_imin"]

# path2_panel_metadata_csv = snakemake.input["panel_metadata"]
# results_fp = snakemake.output['panel_mod_results']

log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #

# Load immune reference dataset
with open(immune_reference_dataset, 'r') as fp:
    imm_ref = pd.read_csv(fp)

# Load immune methylation data - dss_df_30x_probes_betas
with open(meth_results_fp, 'r') as fp:
    dss_df_30x_probes_betas = pd.read_csv(fp)

# Load methatlas deconv output
with open(deconv_fp, "r") as fd:
    deconv_df = pd.read_csv(fd)

deconv_df.rename(
    {"Unnamed: 0": "Immune cells", "sample_beta_vals": "sample % cell"},
    axis=1,
    inplace=True,
)

# Monocytes = deconv_df[deconv_df['Immune cells']=="Monocytes_EPIC"]['sample % cell'].item()
# Bcells = deconv_df[deconv_df['Immune cells']=="B-cells_EPIC"]['sample % cell'].item()
# CD4_Tcells = deconv_df[deconv_df['Immune cells']=="CD4T-cells_EPIC"]['sample % cell'].item()
# NK_cells = deconv_df[deconv_df['Immune cells']=="NK-cells_EPIC"]['sample % cell'].item()
# CD8_Tcells = deconv_df[deconv_df['Immune cells']=="CD8T-cells_EPIC"]['sample % cell'].item()
# Neutrophils = deconv_df[deconv_df['Immune cells']=="Neutrophils_EPIC"]['sample % cell'].item()
# Bladder = deconv_df[deconv_df['Immune cells']=="Bladder"]['sample % cell'].item()

Monocytes = deconv_df[deconv_df['Immune cells']=="CD14"]['sample % cell'].item()
Bcells = deconv_df[deconv_df['Immune cells']=="CD19"]['sample % cell'].item()
CD4_Tcells = deconv_df[deconv_df['Immune cells']=="CD4_Eff"]['sample % cell'].item()
NK_cells = deconv_df[deconv_df['Immune cells']=="CD56"]['sample % cell'].item()
CD8_Tcells = deconv_df[deconv_df['Immune cells']=="CD8"]['sample % cell'].item()
Tregs = deconv_df[deconv_df['Immune cells']=="Treg"]['sample % cell'].item()
Endothelial = deconv_df[deconv_df['Immune cells']=="Endothelial"]['sample % cell'].item()
Eosinophils = deconv_df[deconv_df['Immune cells']=="Eos"]['sample % cell'].item()
Fibroblasts = deconv_df[deconv_df['Immune cells']=="Fibroblast"]['sample % cell'].item()
Neutrophils = deconv_df[deconv_df['Immune cells']=="Neu"]['sample % cell'].item()
Bladder = deconv_df[deconv_df['Immune cells']=="Bladder"]['sample % cell'].item()

# Get ratios of interest
LMR_ratio = get_LMR(deconv_df, log)
NLR_ratio = get_NLR(deconv_df, log)

log.write(f'Proportion of infiltrating Monocytes is {Monocytes}')
log.write(f'Proportion of infiltrating B-cells is {Bcells}')
log.write(f'Proportion of infiltrating CD4+ T-cells is {CD4_Tcells}')
log.write(f'Proportion of infiltrating CD8+ T-cells is {NK_cells}')
log.write(f'Proportion of infiltrating NK cells is {CD8_Tcells}')
log.write(f'Proportion of infiltrating Neutrophils is {Neutrophils}')
log.write(f'Proportion of infiltrating Tregs is {Tregs}')
log.write(f'Proportion of infiltrating Endothelial cells is {Endothelial}')
log.write(f'Proportion of infiltrating Eosinophils is {Eosinophils}')
log.write(f'Proportion of infiltrating Fibroblasts is {Fibroblasts}')
log.write(f'Proportion of Bladder cells is {Bladder}')

results_dict = {
    "LMR_ratio": LMR_ratio,
    "NLR_ratio": NLR_ratio,
    "Mono_inf": Monocytes,
    "Bcells_inf": Bcells,
    "CD4_inf": CD4_Tcells,
    "CD8_inf": CD8_Tcells,
    "NK_inf": NK_cells,
    "Neut_inf": Neutrophils,
    "Bladder_inf": Bladder,
}

with open(results_json_imin, 'w') as fw:
    json.dump(results_dict, fw)
    

log.close()

# Notes:

# Platelet isn't available as they don't have gDNA. Future investigations may reveal the ability to use mDNA signatures for this
# https://www.sciencedirect.com/science/article/pii/S0167527323016960
# Th1 isn't available with this reference
# Dendritic cells aren't available with this reference
# M2 macrophages specifically aren't available with this reference
