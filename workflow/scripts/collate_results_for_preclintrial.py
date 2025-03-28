#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# collate_results_for_preclintrial.py
################################################################################

import pandas as pd
from snakemake.script import snakemake
from shared_functions import VARIANT_TYPE, preclin_stage_panel_result_header, variant_prep

# DEBUG PATHS
    # path2_panel_metadata_csv = f"{PREFIX}/config/testing_panel_metadata1.csv"
    # snv_panel_csv = f"{PREFIX}/results_debug/{SAMPLE}/snv_annotation/{SAMPLE}.snv_results.csv"
    # panel_mod_results = f"{PREFIX}/results/{SAMPLE}/mod_calling/{SAMPLE}.mod_results.csv"
    # immune_results = f"{PREFIX}/results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
    # panel_results = f"{PREFIX}/results_debug/{SAMPLE}/panel_results.csv"


path2_panel_metadata_csv = snakemake.input["panel_metadata"]
snv_results = snakemake.input["snv_results"]
mod_results = snakemake.input["mod_results"]
immune_results = snakemake.input["immune_results"]
panel_results = snakemake.output["panel_results"]


with open(snv_results, 'r') as f1:
    snv_df = pd.read_csv(f1, dtype={"ID": str})
with open(mod_results, 'r') as f1:
    mod_df = pd.read_csv(f1, dtype={"ID": str})
with open(immune_results, 'r') as f1:
    immune_df = pd.read_csv(f1, dtype={"ID": str})

# load and process demographic and clinicopath results
# path2_panel_metadata_csv = snakemake.input["panel_metadata"]

preclin_stage_panel_result_header = ["ID", "Marker name", "Scoring Type", "Biomarker Type", "Result Options", "Result"]
preclin_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)

panel_demog = variant_prep(path2_panel_metadata_csv, variant_type='demographic')
panel_clinpath = variant_prep(path2_panel_metadata_csv, variant_type='clinicopathology')

demclin_df = pd.concat([panel_demog, panel_clinpath])
for i, row in demclin_df.iterrows():
    preclin_panel_df.loc[i] = [row['ID'], row['Gene name'], row['Scoring Type'], row[VARIANT_TYPE], row['Variant'], row['Result']]  # type: ignore  

full_panel_results_df = pd.concat([snv_df, mod_df, immune_df, preclin_panel_df])
full_panel_results_df = full_panel_results_df.sort_values(['ID'])
full_panel_results_df.reset_index(drop = True, inplace = True)

# save preclin panel
full_panel_results_df.to_csv(panel_results, index=False)
