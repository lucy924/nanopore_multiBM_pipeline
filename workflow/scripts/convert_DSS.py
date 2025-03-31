#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# modification_calling.py
################################################################################

import pandas as pd
from snakemake.script import snakemake

# Previous shell code
# awk -v OFS='\t' 'BEGIN{print "chr","pos","N","X"}{print $1,$2,($12+$13),$13}' ${snakemake_input[0]} > ${snakemake_output[0]}
# this is wrong!!!! 12 and 13 should be swapped around!!!


# convert to DSS is getting N and X from bedmethyl
# chrom, start, end, mod base code, score, strand, start pos, end pos, color, N valid cov, fraction modified, N mod, N canonical, N other mod, N delet, N fail, N diff, N nocall

# path2bed = "../ext_bm_pipeline_dev/results/test2/mod_calling/test2.wf_mods.all.bedmethyl.bed"
header = ["chrom", "start", "end", "mod base code", "score", "strand", "start pos", "end pos", "color", "N valid cov", "fraction modified", "N mod", "N canonical", "N other mod", "N delete", "N fail", "N diff", "N nocall"]

with open(snakemake.input[0], 'r') as fb:
    bed_df = pd.read_csv(fb, sep = '\t', names=header)

# bed_df1 = bed_df[bed_df['mod base code'] == 'm']
dss_df = bed_df.copy()[['chrom', 'start', 'N valid cov', 'N mod']]
dss_df.rename({'chrom':'chr', 'start': 'pos', 'N valid cov': "N", 'N mod': "X"}, axis=1, inplace=True)

dss_df.to_csv(snakemake.output[0], sep = '\t', index=False)
