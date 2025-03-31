#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# modification_calling.py
################################################################################


import pandas as pd
from snakemake.script import snakemake
from shared_functions import filter_to_variant_type, add_functional_flanking_regions
from shared_functions import CHROMOSOMES


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
        right_on=["chr", "endpos"],
        how="right",
    )

    merged_df_clean = merged_df[merged_df["N"].notna()].copy()
    merged_df_clean["endpos"] = merged_df_clean["endpos"].astype(int)
    merged_df_clean["N"] = merged_df_clean["N"].astype(int)
    merged_df_clean["X"] = merged_df_clean["X"].astype(int)
    merged_df_clean.drop("start", axis=1, inplace=True)

    # Select the desired columns and rename if necessary
    df = merged_df_clean[["probe", "chr",
                          "endpos", "strand", "N", "X"]]
    return df


def calculate_beta_value(methylated_reads, unmethylated_reads, epsilon=100):
    """
    Calculate the methylation beta value, based on minfi

    Parameters:
    methylated_reads (int): Number of methylated reads (M).
    unmethylated_reads (int): Number of unmethylated reads (U).
    epsilon (float): Constant to avoid division by zero (often set to 100).

    Returns:
    float: The calculated beta value.
    """
    # TODO: check this result to see if it matches minfi's beta values
    beta_value = methylated_reads / (methylated_reads + unmethylated_reads + epsilon)
    return beta_value


def add_beta_vals(df):
    # Derive column 'beta' from columns: 'N', 'X'
    df["beta"] = calculate_beta_value(df["X"], df["N"])
    # Change column type to float64 for column: 'beta'
    df = df.astype({"beta": "float64"})
    return df


def export_for_methatlas(df):
    # Filter rows based on column: 'probe'
    df = df[df['probe'].notna()]
    # Drop columns: 'chr', 'start' and 4 other columns
    df = df.drop(columns=["chr", "endpos", "strand", "N", "X"])
    # Rename column 'beta' to 'LPPP_Tumour_beta_vals'
    df = df.rename(columns={"beta": "sample_beta_vals"})

    # fp = f"{sample_name}.methatlas.csv"
    fp = snakemake.output['epic_probe_results']
    df.to_csv(fp, index=False)

    return


def calculate_prom_start(row):
    # Custom function to calculate prom_start
    if row["Is this record the whole gene?"] == "Yes":
        if row["strand"] == "+":
            return row["start pos"] - 2000
        else:
            return row["end pos"]
    else:
        return row["start pos"] - 5


def calculate_prom_end(row):
    # Custom function to calculate prom_end
    if row["Is this record the whole gene?"] == "Yes":
        if row["strand"] == "+":
            return row["start pos"]
        else:
            return row["end pos"] + 2000
    else:
        return row["end pos"] + 5
        
        
def add_prom_start_end(panel_metadata_df_meth):
    # Custom function to calculate prom_start
    # note that if the entry is not the whole gene it makes it flank the region to ensure the row and result gets included later

    # Apply the custom function to create the "prom_start" and "prom_end" column
    panel_metadata_df_meth["prom_start"] = panel_metadata_df_meth.apply(
        calculate_prom_start, axis=1
    )
    panel_metadata_df_meth["prom_end"] = panel_metadata_df_meth.apply(
        calculate_prom_end, axis=1
    )

    return panel_metadata_df_meth


def calculate_down_start(row):
    if row["Is this record the whole gene?"] == "Yes":
        if row["strand"] == "+":
            return row["end pos"]
        else:
            return row["start pos"] - 1000
    else:
        return row["start pos"] - 5


def calculate_down_end(row):
    # Custom function to calculate prom_end
    if row["Is this record the whole gene?"] == "Yes":
        if row["strand"] == "+":
            return row["end pos"] + 1000
        else:
            return row["start pos"]
    else:
        return row["end pos"] + 5
        
        
def add_downstream_start_end(panel_metadata_df_meth):
    # note that if the entry is not the whole gene it makes it flank the region to ensure the row and result gets included later

    # Apply the custom function to create the "prom_start" and "prom_end" column
    panel_metadata_df_meth["down_start"] = panel_metadata_df_meth.apply(
        calculate_down_start, axis=1
    )
    panel_metadata_df_meth["down_end"] = panel_metadata_df_meth.apply(
        calculate_down_end, axis=1
    )

    return panel_metadata_df_meth


def meth_reduce_metadata_df(df):
    """Drop columns of excess info - for methylation only"""
    df = df.drop(
        columns=[
            "Is this record the whole gene?",
            "BCG response characteristic",
            "BCG non response characteristic",
            "Is variant in coding region?",
            "SNP ID",
            "Result Options",
            "Illumina EPIC ID",
            "DNA methylation",
            "Notes",
            "WARNINGS",
            "Ref",
            "length",
        ]
    )
    return df


def sort_and_prep_dfs(variants_df, dss_betas_df):
    # set chromosome column types
    variants_df.chrom = variants_df.chrom.astype("category")
    variants_df.chrom = variants_df.chrom.cat.set_categories(CHROMOSOMES)
    dss_betas_df.chr = dss_betas_df.chr.astype("category")
    dss_betas_df.chr = dss_betas_df.chr.cat.set_categories(CHROMOSOMES)

    print(variants_df)
    # convert locations to ints
    variants_df["prom_start"] = variants_df["prom_start"].astype("int64")
    variants_df["prom_end"] = variants_df["prom_end"].astype("int64")
    variants_df["down_start"] = variants_df["down_start"].astype("int64")
    variants_df["down_end"] = variants_df["down_end"].astype("int64")
    variants_df["start pos"] = (
        variants_df["start pos"].str.replace(",", "").astype("int64")
    )
    variants_df["end pos"] = (
        variants_df["end pos"].str.replace(",", "").astype("int64")
    )
    # dss_betas_df["start"] = dss_betas_df["start"].astype("int64")
    dss_betas_df["endpos"] = dss_betas_df["endpos"].astype("int64")

    # Sort by columns: 'chrom' (ascending), 'start pos' (ascending), 'end pos' (ascending)
    variants_df = variants_df.sort_values(["chrom", "start pos", "end pos"])
    dss_betas_df = dss_betas_df.sort_values(["chr", "endpos"])
    return variants_df, dss_betas_df


def merge_data_to_list_dicts(results_df, panel_data):
    results_as_list_dicts = list()
    for i in results_df.groupby(by="chr", observed=False):
        chrom = i[0]
        print(f"Processing chromosome: {chrom}")
        if chrom not in set(panel_data["chrom"]):
            print("    Chromosome not in panel regions of interest; skipping.")
            continue

        # Get chromosome group from variants
        df = panel_data.copy()
        print(df.head())
        chr_entries = df.loc[df["chrom"] == i[0]].copy()
        print(chr_entries.head())
        chr_entries["min"] = chr_entries[
            ["prom_start", "prom_end", "down_start", "down_end"]
        ].min(axis=1)
        chr_entries["max"] = chr_entries[
            ["prom_start", "prom_end", "down_start", "down_end"]
        ].max(axis=1)

        # check if endpos is in any of the metadata entries
        # if so, add all data to results
        for i, row_a in i[1].iterrows():
            for j, row_b in chr_entries.iterrows():
                if row_b["min"] <= row_a["endpos"] <= row_b["max"]:
                    merged_dict = {**row_a.to_dict(), **row_b.to_dict()}
                    cpg_loc = merged_dict["start"]
                    merged_dict["CpG location"] = cpg_loc
                    keys = ["chr", "endpos", "start", "min", "max"]
                    list(map(merged_dict.pop, keys))
                    merged_dict["CpG region"] = ""
                    if (
                        merged_dict["start pos"] <= cpg_loc <= merged_dict["end pos"]
                    ) or (
                        merged_dict["start pos"] >= cpg_loc >= merged_dict["end pos"]
                    ):
                        merged_dict["CpG region"] = "intragenic"
                    elif (
                        merged_dict["prom_start"] <= cpg_loc <= merged_dict["prom_end"]
                    ) or (
                        merged_dict["prom_start"] >= cpg_loc >= merged_dict["prom_end"]
                    ):
                        merged_dict["CpG region"] = "promoter"
                    elif (
                        merged_dict["down_start"] <= cpg_loc <= merged_dict["down_end"]
                    ) or (
                        merged_dict["down_start"] >= cpg_loc >= merged_dict["down_end"]
                    ):
                        merged_dict["CpG region"] = "downstream"
                    results_as_list_dicts.append(merged_dict)

    return results_as_list_dicts


def make_results_df(results_as_list_dicts):
    results_df = pd.DataFrame(
        results_as_list_dicts,
        columns=[
            "ID",
            "Biomarker name",
            "chrom",
            "start pos",
            "end pos",
            "strand",
            "probe",
            "beta",
            "CpG location",
            "CpG region",
            "N",
            "X",
            "prom_start",
            "prom_end",
            "down_start",
            "down_end",
        ],
    )

    # Change column type to int64 for columns: 'start pos', 'end pos' and 5 other columns
    results_df = results_df.astype(
        {
            "start pos": "int64",
            "end pos": "int64",
            "CpG location": "int64",
            "prom_start": "int64",
            "prom_end": "int64",
            "down_start": "int64",
            "down_end": "int64",
        }
    )

    # Sort by columns: 'id_num' (ascending), 'CpG location' (ascending)
    results_df["id_num"] = results_df["ID"].astype("int64")
    results_df = results_df.sort_values(["id_num", "CpG location"])
    results_df.drop("id_num", axis=1, inplace=True)

    results_df.reset_index(inplace=True, drop=True)

    return results_df


def get_results_df(results_df, panel_data):
    """these steps to merge dfs are broken up because of the long time it takes to go through the dss results"""
    results_as_list_dicts = merge_data_to_list_dicts(
        results_df, panel_data
    )
    results_df = make_results_df(results_as_list_dicts)
    return results_df

# ================ Process probes for immune infiltrate ================ #
# ------------------------------------------------ #
# import EPIC genomic locations
# EPIClocs_fp = "../IlluminaEPIC_genomic_locations_hg38.csv"
EPIClocs_fp = snakemake.input['IlluminaEPIC_genomic_locations_hg38']
with open(EPIClocs_fp, "r") as fp:
    EPIClocs_df = pd.read_csv(fp, sep=",", index_col=0)
    
# ------------------------------------------------ #
# Prep wf output file
# dss_fp = os.path.join(data_dir, sample_name  + ".wf_mods.all.dss_format.tsv")
dss_fp = snakemake.input["dss_file"]
with open(dss_fp, "r") as fp:
    dss_df = pd.read_csv(fp, sep="\t")
    
dss_df_30x = filter_dss_to_x_depth(dss_df, depth=1)  # default 30x
#TODO: remove depth filter

# add illumina probes where applicable, keep all dss data though
dss_df_30x_probes = add_illumina_probes(dss_df_30x, EPIClocs_df)
dss_df_30x_probes_betas = add_beta_vals(dss_df_30x_probes)


# ------------------------------------------------ #
# Export data for immune infiltrate

export_for_methatlas(dss_df_30x_probes_betas)

# ================ Process panel methylation ================ #

# Load in metadata files
# hans_prefix = "/home/dejlu879/"
# path2_panel_metadata_csv = hans_prefix + "ProjectProtocol/20231211-prep_for_ESR/20231211-metadata.csv"
panel_csv_fp = snakemake.input["panel_metadata"]

# with open(path2_panel_metadata_csv, "r") as fpc:
#     panel_metadata_df = pd.read_csv(fpc, dtype={"ID": str})

# Load panel
with open(panel_csv_fp, 'r') as fp:
    panel_csv = pd.read_csv(fp, dtype={'ID': str}, thousands = ',')

# ------------------------------------------------ #
# Locate genomic coordinates of key points (from metadata file)
panel_metadata_df = filter_to_variant_type(panel_csv, variant_type = 'mod')

# panel_meth_flank_df = add_prom_start_end(panel_metadata_df_meth.copy())
# panel_meth_flank_df = add_downstream_start_end(
#     panel_meth_flank_df.copy())

# Add flanking regions
panel_csv = add_prom_start_end(panel_csv.copy())
panel_csv = add_downstream_start_end(panel_csv.copy())
# Note: don't replace these with add_functional_flanking_regions because the distinctions are neccessary for categorisation of the variant results

# ------------------------------------------------ #
# Compare methylation discovered at those coords to what is expected/useful (from metadata file)
# meth_threshold = 0.8 (? discuss with others)  
# `pos + 1` should map correctly to the methylated loci

panel_metadata_reduced_df = meth_reduce_metadata_df(
    panel_metadata_df)

variants_meth_reduced_sorted_df, dss_df_30x_probes_betas_sorted = (
    sort_and_prep_dfs(
        panel_metadata_reduced_df, dss_df_30x_probes_betas
    )
)

# ------------------------------------------------ #
# Create a new dataframe containing all the cpgs in each desired location, getting all info from variants into each entry in results
results_df = get_results_df(
    dss_df_30x_probes_betas_sorted, variants_meth_reduced_sorted_df
)

results_df.to_csv(snakemake.output['panel_mod_results'])

