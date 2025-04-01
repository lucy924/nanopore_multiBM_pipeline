#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# generate_report_template.py
# see ProjectProtocolDev -> Sort out scoring code #2
# Starting from https://stackoverflow.com/a/63622717
################################################################################

from jinja2 import Template
import codecs
# from snakemake.script import snakemake
import json
import pandas as pd
from shared_functions import BIOMARKER_TYPE


import os
if os.getenv("SNAKEMAKE_DEBUG"):
    class FakeSnakemake:
        SAMPLE = "test4predict"
        PROJECT = "BCG_on_NMIBC"
        PREFIX = "/external/analyses/lucy/nanopore_multiBM_pipeline"
        input = {
            'panel_metadata': f"{PREFIX}/config/testing_panel_metadata_w_scores.csv",
            'report_template': f"{PREFIX}/resources/template.md",
            'scoring_results': f"{PREFIX}/results_debug/{SAMPLE}/{SAMPLE}.scores.csv"
            }
        output = {
            'report_md': f"{PREFIX}/results_debug/{SAMPLE}/{SAMPLE}.report.md"
            }
        params = {
            'sample_name': SAMPLE,
            'panel_name': PROJECT
            }
        log = [f"{PREFIX}/results_debug/{SAMPLE}.report_template.log"]

    snakemake = FakeSnakemake()

# ------------------------------------------------ #
# # get snakemake variables
sample_name = snakemake.params['sample_name']
panel_name = snakemake.params['panel_name'].replace('_', ' ')

path2_panel_metadata_csv = snakemake.input["panel_metadata"]
path2_template = snakemake.input["report_template"]
scores_csv = snakemake.input["scoring_results"]
results_fp = snakemake.output['report_md']

log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #

# create a dict with all data that will populate the template
    
# Get results from csv into dict
with open(scores_csv, 'r') as fp:
    results_df = pd.read_csv(fp, index_col=0)
with open(path2_panel_metadata_csv, 'r') as fp:
    panel_df = pd.read_csv(fp, index_col=0, dtype={'ID': str})
    panel_df.rename(columns = {BIOMARKER_TYPE: 'Biomarker Type'}, inplace=True)
    
header = ["ID", "Biomarker name", "Biomarker Type", "Result Options", "Result", "Score"]
panel_cols = ["Biomarker name", "Biomarker Type", "Result Options"]

# Merge using index
all_results_df = results_df.merge(panel_df[panel_cols], left_index=True, right_index=True, how='left', )
all_results_df.reset_index(names='ID', inplace=True)
all_results_df = all_results_df[header]
    
# all_results_dict = all_results_df.to_dict()
all_results_df = all_results_df.map(lambda x: x.replace('|', '\\|') if isinstance(x, str) else x)

results = {
    "sample_name": sample_name,
    "panel_name": panel_name,
    "url": "https://github.com/lucy924/nanopore_multiBM_pipeline",
}

# results["all_results"] = list_of_lists = [
#     [key, value] for key, value in all_results_dict['Score'].items()
# ]

all_data_lists = list()
for idx, entry in all_results_df.iterrows():
    # entry_list = [idx]
    # entry_list.extend(list(entry))
    all_data_lists.append(list(entry))
    
results["all_results"] = all_data_lists

# get word result from normalised score
normalised_score = all_results_df['Score'].iloc[-1]
if normalised_score < 25:
    results["word_result"] = 'Highly Likely'
elif normalised_score < 50:
    results["word_result"] = 'Somewhat Likely'
elif normalised_score < 75:
    results["word_result"] = 'Somewhat Unlikely'
else:
    results["word_result"] = 'Highly Unlikely'
    

# render the template
# path2template = (
    # "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/resources/template.md"
# )

with open(path2_template, "r") as file:
    template = Template(file.read(), trim_blocks=True)
rendered_file = template.render(repo=results)

# output the file
output_file = codecs.open(results_fp, "w", "utf-8")
output_file.write(rendered_file)
output_file.close()

log.close()
