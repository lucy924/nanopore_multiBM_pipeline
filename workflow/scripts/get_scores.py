#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# get_scores.py
################################################################################


import os
import json
import numpy as np
from numpy import float64
import pandas as pd
from pprint import pprint
from snakemake.script import snakemake


def match_genotype(allele1, allele2, vals, hrsd_vals, scored):
    # THERE MAY BE MULTIPLE OPTIONS
    # split options and zip them with associated HR/SD
    vals_options = vals.split('/')
    hrsd_options = hrsd_vals.split('/')
    vals_with_hrsd = zip(vals_options, hrsd_options)
    for vals_genotype, hrsd in vals_with_hrsd:
        panel_allele1, panel_allele2 = vals_genotype.split('|')
        
        # section that matches allele options
        if allele1 == panel_allele1 and allele2 == panel_allele2:
            score = float64(hrsd)
            scored = True
            return score, scored
        elif allele1 == panel_allele2 and allele2 == panel_allele1:
            score = float64(hrsd)
            scored = True
            return score, scored

# ------------------------------------------------ #
# get snakemake variables
path2_panel_metadata_csv = snakemake.input["panel_metadata"]
log = open(snakemake.log[0], 'w')
# ------------------------------------------------ #

with open(path2_panel_metadata_csv, "r") as fpc:
    panel_metadata_df = pd.read_csv(fpc, dtype={"ID": str})

panel_metadata_df = panel_metadata_df[
    (panel_metadata_df["HR/SD Response"] == variant_type)
]
    



log.close()
