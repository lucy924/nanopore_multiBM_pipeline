#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# get_patient_results.py
################################################################################


import os
import json
import numpy as np
import pandas as pd
from snakemake.script import snakemake

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
