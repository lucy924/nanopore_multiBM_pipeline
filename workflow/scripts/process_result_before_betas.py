#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# process_result_before_betas.py
################################################################################


import pandas as pd
from snakemake.script import snakemake  # type: ignore


def filter_dss_to_x_depth(dss_df, depth = 30):
    # Filter rows based on column: 'N'
    dss_df = dss_df[dss_df["N"] >= depth]
    # Created column 'endpos' from formula
    dss_df["endpos"] = dss_df["pos"] + 1
    # dss_df['b4pos'] = dss_df['pos']-1
    return dss_df



def add_illumina_probes(dss_df, EPIClocs_df):
    # Merge df2 and df1 based on seqnames and start
    merged_df = pd.merge(
        EPIClocs_df,
        dss_df,
        left_on=["seqnames", "start"],
        right_on=["chr", "pos"],
        how="right",
    )

    merged_df_clean = merged_df[merged_df["N"].notna()].copy()
    merged_df_clean["pos"] = merged_df_clean["pos"].astype(int)
    merged_df_clean["N"] = merged_df_clean["N"].astype(int)
    merged_df_clean["X"] = merged_df_clean["X"].astype(int)
    merged_df_clean.drop("start", axis=1, inplace=True)

    # Select the desired columns and rename if necessary
    df = merged_df_clean[["probe", "chr",
                          "pos", "strand", "N", "X"]]
    return df


def export_for_getting_betas(df):
    # Filter rows based on column: 'probe'
    # df = df[df['probe'].notna()]
    # Drop columns: 'chr', 'start' and 4 other columns
    # df = df.drop(columns=["chr", "start", "endpos", "strand", "N", "X"])
    # df = df.drop(columns=["chr", "endpos", "strand", "N", "X"])
    # Rename column 'beta' to 'LPPP_Tumour_beta_vals'
    # df = df.rename(columns={"beta": "sample_beta_vals"})

    # fp = f"{sample_name}.methatlas.csv"
    fp = snakemake.output['epic_probe_results']
    df.to_csv(fp, index=False)

    return

# ------------------------------------------------ #
# get snakemake variables
EPIClocs_fp = snakemake.input['IlluminaEPIC_genomic_locations_hg38']
# EPIClocs_fp = "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/resources/IlluminaEPIC_genomic_locations_hg38.csv"
dss_fp = snakemake.input["dss_file"]
# dss_fp = "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/results/test2/mod_calling/test2.wf_mods.all.dss_format.tsv"
# path2_panel_metadata_csv = snakemake.input["panel_metadata"]
# results_fp = snakemake.output['panel_mod_results']

log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #
# import EPIC genomic locations
# EPIClocs_fp = "../IlluminaEPIC_genomic_locations_hg38.csv"
with open(EPIClocs_fp, "r") as fp:
    EPIClocs_df = pd.read_csv(fp, sep=",")
    
# ------------------------------------------------ #
# Prep wf output file
# dss_fp = os.path.join(data_dir, sample_name  + ".wf_mods.all.dss_format.tsv")
with open(dss_fp, "r") as fp:
    dss_df = pd.read_csv(fp, sep="\t")
    
dss_df_30x = filter_dss_to_x_depth(dss_df.copy(), 1)  #TODO: default 30x, use 1x for trial

# add illumina probes where applicable, keep all dss data though
dss_df_30x_probes = add_illumina_probes(dss_df_30x, EPIClocs_df)


# ------------------------------------------------ #
# Export data for immune infiltrate

fp = snakemake.output['pre_beta']
dss_df_30x_probes.to_csv(fp, index=False)
