#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# calculate_coverage.py
################################################################################

import os

if os.getenv("SNAKEMAKE_DEBUG"):
    class FakeSnakemake:
        SAMPLE = "findorig_adaptiverefcov"
        PROJECT = "sequenced_adaptiveref"
        BUFFER = 2000
        PREFIX = "/external/analyses/lucy/nanopore_multiBM_pipeline"
        input = {
            'minknow_target_bed': f"results/{PROJECT}/minknow_input_supp/targets.minknow.{BUFFER}.bed"
            }
        output = {
            'final_bed_name': f"results/{PROJECT}/minknow_input/targets.buffed.bed"
            }
        params = {
            'min_cov': 0.5,
            'max_cov': 2.0,
            'buffer': BUFFER
        }
        log = [f"{PREFIX}/results_debug/{SAMPLE}.snv_annotation.log"]

    snakemake = FakeSnakemake()

import pandas as pd
import json
# from snakemake.script import snakemake  # type: ignore
from shared_functions import HG_LENGTH, CHROMOSOMES, CHROM_LENGTHS

def calc_coverage(variants_bed, log):
    
    chr_coverage = dict()
    total_length = 0

    for rowid, entry in variants_bed.iterrows():
        chr_name = entry['#chrom']
        chrStart = entry['chromStart']
        chrEnd = entry['chromEnd']
        if chr_name not in chr_coverage:
            chr_coverage[chr_name] = 0
        chr_coverage[chr_name] += int(chrEnd) - int(chrStart)
        total_length += int(chrEnd) - int(chrStart)

    for chr in CHROMOSOMES:
        if chr not in chr_coverage:
            chr_coverage[chr] = 0
        log.write(f'total length in {chr}: {chr_coverage[chr]}\n')
        log.write(f'percent coverage: {(chr_coverage[chr]/CHROM_LENGTHS[chr])*100}\n\n')
            
    print(f'total length of areas in bed file {total_length:,} bp')
    log.write(f'total length of areas in bed file {total_length:,} bp\n')
    perc_HG = round((total_length/HG_LENGTH)*100, 2)
    print(f'percent of HG: {perc_HG}%')
    log.write(f'percent of HG: {perc_HG}%\n')
    
    return perc_HG

log = open(snakemake.log[0], 'w')
input_bed = snakemake.input['minknow_target_bed']
bed_df = pd.read_csv(input_bed, sep = '\t', names = ['#chrom', 'chromStart', 'chromEnd', 'name'])
perc_HG = calc_coverage(bed_df, log)

# Check paramters
min_cov = snakemake.params['min_cov']
max_cov = snakemake.params['max_cov']
buffer = snakemake.params['buffer']

log.write('--------------------------------\n')
output_bed_fp = snakemake.output['final_bed_name']

if min_cov < perc_HG < max_cov:
    log.write(f"Yay we found it! Final coverage: {perc_HG}%\n")
    print(f"Yay we found it! Final coverage: {perc_HG}%")
    os.system(f'cp {input_bed} {output_bed_fp}')
    # "Criteria met with value {config[start_value]}" > {output}
else:
    log.write(
        "Criteria not met. Please adjust 'flanking_bp' in config/config.yaml and rerun.\n" +
        f"Selected flanking_bp ({buffer}) gave percentage coverage: {perc_HG}%\n" +
        f"Goal: between {min_cov} and {max_cov}"
        )
    print(
        "Criteria not met. Please adjust 'flanking_bp' in config/config.yaml and rerun.\n",
        f"Selected flanking_bp ({buffer}) gave percentage coverage: {perc_HG}%\n",
        f"Goal: between {min_cov} and {max_cov}"
        )
    with open(output_bed_fp, 'w') as fw:
        fw.write('Criteria not met.')

log.close()
