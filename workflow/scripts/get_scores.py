#!/usr/bin/env python
# -*- mode: python; python-indent: 4; tab-width: 4; indent-tabs-mode: nil -*-
################################################################################
# LUCY PICARD
# ---------------------------------- #
# get_scores.py
# see ProjectProtocolDev -> Sort out scoring code #2
################################################################################

import os
if os.getenv("SNAKEMAKE_DEBUG"):
    class FakeSnakemake:
        SAMPLE = "test4predict"
        PREFIX = "/external/analyses/lucy/nanopore_multiBM_pipeline"
        input = {
            'panel_metadata': f"{PREFIX}/config/testing_panel_metadata_w_scores.csv",
            'snv_results': f"{PREFIX}/results/{SAMPLE}/snv_annotation/{SAMPLE}.snv_results.csv",
            'mod_results': f"{PREFIX}/results/{SAMPLE}/mod_calling/{SAMPLE}.mod_results.csv",
            'immune_results': f"{PREFIX}/results/{SAMPLE}/immune_infiltrate/{SAMPLE}.immune_panel_results.csv"
            }
        output = {
            'scores': f"{PREFIX}/results_debug/{SAMPLE}/{SAMPLE}.scores.csv"
            }
        log = [f"{PREFIX}/results_debug/{SAMPLE}.get_scores.log"]

    snakemake = FakeSnakemake()


import os
import json
import numpy as np
from numpy import float64
import pandas as pd
# from snakemake.script import snakemake
from shared_functions import BIOMARKER_TYPE


# functions

# VARIANT_TYPE = "Variant Type (snv, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, demographic, clinicopathology)"

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


def calculate_score(control_val, result_val, std_dev, HRSD_val, response):
    
    if response:
        # TODO: this needs an expert pair of eyes on it
        # Calculations between 0 - 1 - [2 and up] for hazard ratio means I don't know how best to calculate this for going down om score. Probably depends on actual results from the classifier too.
        # find the deviation away from 1, don't use the actual value
        HRSD_val = (1 - HRSD_val)
    else:  # non-response
        HRSD_val = (HRSD_val - 1)
        
    # Number of std_dev's between control and result"""
    num_std_devs = abs(control_val - result_val) / std_dev
    # What is the final score given this many standard deviations? in %
    score = (HRSD_val * num_std_devs)  #*100  
    # if nonresponse, add the 1 back
    if not response:
        score = (score + 1)
    
    return score
        
        
def calculate_continuous_score(result, SD, control_val, HRSD, group):
    # result = 42
    if group == 'high':
        HRSD_dif = HRSD - 1
    elif group == 'low':
        HRSD_dif = 1 - HRSD
    else:
        raise ValueError("group needs to be specified.")
    
    score = calculate_score(control_val, result, SD, HRSD_dif, response = False)
    # note that response is only false here because I've already taken the difference into account
    
    # # First how many std devs are there between the mean and the result?
    # dif_btwn_result_and_mean = result - control_val
    # num_stddevs = dif_btwn_result_and_mean / SD
    # # now what is the percentage score?
    # score = (HRSD_dif * num_stddevs)*100
    
    return score


def get_cont_score_when_both_R_NR(panel_entry, result_entry):
    """Scoring when there is continuous data for both response and non-response, and there needs to be a cap in the middle so the two sides don't overlap"""
    # TODO: check this scoring makes sense, particularly with high/low/resp/nonresp thresholds
    
    # Get HR/SD for high and low groups
    high_HRSD = float64(panel_entry['HR/SD Non-Response'])
    low_HRSD = float64(panel_entry['HR/SD Response'])
    # high_HRSD = 1.17
    # low_HRSD = 0.83
    
    result = result_entry['Result']  # 0.524186
    control_group_value = float64(panel_entry['Control Group or Baseline'])
    SD_value_resp = float64(panel_entry['Response StdDev'])
    SD_value_nonresp = float64(panel_entry['Non-Response StdDev'])
    high_threshold = control_group_value + SD_value_nonresp
    low_threshold = control_group_value - SD_value_resp
    
    # Find out if the result belongs to Response or Non-Response
    high_group = False
    low_group = False
    if result >= high_threshold:
        print(f'result {result} is in the high group')
        high_group = True
    elif result <= low_threshold:
        print(f'result {result} is in the low group')
        low_group = True
    else:
        print(f"result {result} is not in any group, don't count it")
        return None
    
    # get appropriate HR/SD value
    if high_group:
        score = calculate_continuous_score(result, SD = SD_value_nonresp, control_val = control_group_value, HRSD = high_HRSD, group = 'high')
    elif low_group:
        score = calculate_continuous_score(result, SD = SD_value_resp, control_val = control_group_value, HRSD = low_HRSD, group = 'low')
    
    print(score)
    return score


def is_response_or_nonresponse(panel_entry):
    response=False
    nonresponse=False
    if panel_entry['Is there response data?'] == 'Yes':
        response = True
    if panel_entry['Is there non-response data?'] == 'Yes':
        nonresponse = True
    return response, nonresponse


def get_scoring_response_from_panel(panel_entry):

    response, nonresponse = is_response_or_nonresponse(panel_entry)
    # high_HRSD = 1.17
    # low_HRSD = 0.83
    
    if response and nonresponse:
        # return a list for the value and the SD, response then nonresponse
        # score = get_cont_score_when_both_R_NR(panel_entry, result_entry)
        # return score
        HRSD_resp = str(panel_entry['HR/SD Response'])
        SD_value_resp = str(panel_entry['Response StdDev'])
        HRSD_nonresp = str(panel_entry['HR/SD Non-Response'])
        SD_value_nonresp = str(panel_entry['Non-Response StdDev'])
        return [HRSD_resp, HRSD_nonresp], [SD_value_resp, SD_value_nonresp], response, nonresponse
            
    if response:
        HRSD_resp = str(panel_entry['HR/SD Response'])
        SD_value = str(panel_entry['Response StdDev'])
    elif nonresponse:
        HRSD_resp = str(panel_entry['HR/SD Non-Response'])
        SD_value = str(panel_entry['Non-Response StdDev'])
    else:
        raise ValueError("uh this shouldn't happen, got no response data apparently")
    
    return HRSD_resp, SD_value, response, nonresponse


# def score_continuous_data(panel_entry, result_entry, HRSD_resp, SD_value):
    #     """"""
        
    #     response, nonresponse = is_response_or_nonresponse(panel_entry)
            
    #     result = float64(result_entry['Result'])  # 0.524186
    #     control_group_value = float64(panel_entry['Control Group or Baseline'])
    #     # control_group_value = 0.65
    #     # SD_value = 0.1
        
    #     # ------ get appropriate HR/SD value ------ #
    #     # is result higher or lower than control group?
    #     higher = False
    #     lower = False
    #     if result == control_group_value:
    #         score = 1.0
    #         return score
    #     else:
    #         score = calculate_score(control_group_value, result, SD_value, HRSD_resp)
    #     # elif result > control_group_value:
    #     #     higher = True
    #     # elif result < control_group_value:
    #     #     lower = True
        
    #     return score


def get_min_max_scores(panel_entry):
    response, nonresponse = is_response_or_nonresponse(panel_entry)
    scoring_type = panel_entry['Scoring Type']
    result_options = panel_entry['Result Options']
    
    if scoring_type == 'genotypic':
        # This structure covers the possibility of having both response and nonresponse data
        min_score = None
        max_score = None
        if response:
            HR_scores = panel_entry['HR/SD Response'].split('/')
            HR_scores = np.array(HR_scores, dtype=np.float64)  # Convert to np.float64
            min_score = min(HR_scores)
        if nonresponse:
            HR_scores = panel_entry['HR/SD Non-Response'].split('/')
            HR_scores = np.array(HR_scores, dtype=np.float64)  # Convert to np.float64
            max_score = min(HR_scores)
        if not min_score:
            min_score = np.float64(1.0)
        if not max_score:
            max_score = np.float64(1.0)    
    
    elif scoring_type == 'continuous':
        min_result, max_result = np.array(result_options.split('-'), dtype=np.float64)  # make it into a numbers list
        if response and nonresponse:
            min_score = calculate_score(
                control_val=float64(panel_entry['Control Group or Baseline']),
                result_val=min_result,
                std_dev=float64(panel_entry['Response StdDev']),
                HRSD_val=float64(panel_entry['HR/SD Response']),
                response = True
                )
            max_score = calculate_score(
                control_val=float64(panel_entry['Control Group or Baseline']),
                result_val=max_result,
                std_dev=float64(panel_entry['Non-Response StdDev']),
                HRSD_val=float64(panel_entry['HR/SD Non-Response']),
                response = False
                )
        elif nonresponse:
            # find the max possible score
            min_score = np.float64(1.0)
            max_score = calculate_score(
                control_val=float64(panel_entry['Control Group or Baseline']),
                result_val=max_result,
                std_dev=float64(panel_entry['Non-Response StdDev']),
                HRSD_val=float64(panel_entry['HR/SD Non-Response']),
                response = False
                )
        elif response:
            max_score = np.float64(1.0)
            min_score = calculate_score(
                control_val=float64(panel_entry['Control Group or Baseline']),
                result_val=min_result,
                std_dev=float64(panel_entry['Response StdDev']),
                HRSD_val=float64(panel_entry['HR/SD Response']),
                response = True
                )
        else:
            raise ValueError("Shouldn't get here")
        
    elif scoring_type == 'binary':
        # note that binary does not have a std dev entry
        # This structure covers the possibility of having both response and nonresponse data
        min_score = None
        max_score = None
        if response:
            HR_score = float64(panel_entry['HR/SD Response'])
            min_score = HR_score
        if nonresponse:
            HR_score = float64(panel_entry['HR/SD Non-Response'])
            max_score = HR_score
        if not min_score:
            min_score = np.float64(1.0)
        if not max_score:
            max_score = np.float64(1.0)    
    else:
        raise ValueError(f'Need to figure out this scoring type in min-maxing: {scoring_type}')
        
    
    return min_score, max_score



# ------------------------------------------------ #
# get snakemake variables
path2_panel_metadata_csv = snakemake.input["panel_metadata"]
snv_results = snakemake.input["snv_results"]
mod_results = snakemake.input["mod_results"]
immune_results = snakemake.input["immune_results"]
scores_csv = snakemake.output["scores"]

log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #

with open(path2_panel_metadata_csv, 'r') as f1:
    # filter panel to stuff that has results
    panel_df = pd.read_csv(f1, dtype={"ID": str})
    panel_df = panel_df[panel_df['Is there either data?'] == 'Yes']
with open(snv_results, 'r') as f1:
    snv_df = pd.read_csv(f1, dtype={"ID": str})
with open(mod_results, 'r') as f1:
    mod_df = pd.read_csv(f1, dtype={"ID": str})
with open(immune_results, 'r') as f1:
    immune_df = pd.read_csv(f1, dtype={"ID": str})


panel_list = panel_df.to_dict(orient='records')
results_snv_list = snv_df.set_index('ID', drop = False).to_dict(orient='index')
results_mod_list = mod_df.set_index('ID', drop = False).to_dict(orient='index')
results_immune_list = immune_df.set_index('ID', drop = False).to_dict(orient='index')
score_results_dict = dict()
min_scores_list = list()
max_scores_list = list()

for panel_entry in panel_list:
    panel_id = panel_entry['ID']
    log.write('-----------------------------')
    log.write(f"Panel ID: {panel_id}")
    # plog.write(entry)
    if panel_id in results_snv_list.keys():
        result_entry = results_snv_list[panel_id]
    elif panel_id in results_mod_list.keys():
        result_entry = results_mod_list[panel_id]
    elif panel_id in results_immune_list.keys():
        result_entry = results_immune_list[panel_id]
    else:
        log.write(f'No data for {panel_id}')
        continue
    
    # log.write(result_entry)
    # break
    
    scoring_type = panel_entry['Scoring Type']
    log.write(scoring_type)
    
    if panel_entry[BIOMARKER_TYPE] == 'immune_inf':
        result_entry['Result'] = result_entry['Result']/100  
        # immune infiltrate is recorded in normal percentage values so messes up scoring
    
    scored = False
    
    if scoring_type == "continuous":
        HRSD_resp, SD_value, response, nonresponse = get_scoring_response_from_panel(panel_entry)
        if response and nonresponse:
            score = get_cont_score_when_both_R_NR(panel_entry, result_entry)
            scored = True
        else:
            control_val = float64(panel_entry['Control Group or Baseline'])
            result = float64(result_entry['Result'])
            if result == control_val:
                score = 1.0
            else:
                if response:
                    score = calculate_score(control_val, result, float64(SD_value), float64(HRSD_resp), response=True)  # type: ignore
                elif nonresponse:
                    score = calculate_score(control_val, result, float64(SD_value), float64(HRSD_resp), response=False)  # type: ignore
                else:
                    raise ValueError("Shouldn't get here")
            # TODO: check
            scored = True
        if scored:
            score_results_dict[panel_id] = (score, result, response, nonresponse)
            if response:
                log.write(f'Response score: {score}')
            if nonresponse:
                log.write(f'NonResponse score: {score}')
        else:
            log.write(f'Entry {panel_id} not scored for some reason')
            continue
        
    elif scoring_type == "genotypic":
        # Get result alleles
        result_allele1, result_allele2 = result_entry['Result'].split('|')
        # log.write(result_allele1)
        # log.write(result_allele2)
        
        # ----- Look for scoring in control data ----- #
        control_vals = panel_entry['Control Group or Baseline']
        if pd.isna(control_vals) or control_vals == '':
            raise ValueError('"Control Group or Baseline" must have an entry.')
        control_vals_list = control_vals.split('/')
        for control_val in control_vals_list:
            control_allele1, control_allele2 = control_val.split('|')
            if result_allele1 == control_allele1 and result_allele2 == control_allele2:
                score = 1.0
                scored = True
            elif result_allele1 == control_allele2 and result_allele2 == control_allele1:
                score = 1.0
                scored = True
        
            if scored:
                log.write(f'score found for {panel_id}')
                break
        
        if scored:
            score_results_dict[panel_id] = (score, result, response, nonresponse)
            if response:
                log.write(f'Response score: {score}')
            if nonresponse:
                log.write(f'NonResponse score: {score}')
            continue
        
        # ----- Look for scoring in response data ----- #
        HRSD_resp, SD_value, response, nonresponse = get_scoring_response_from_panel(panel_entry)
        
        # Get result alleles
        result_allele1, result_allele2 = result_entry['Result'].split('|')
        result = result_entry['Result']
        
        # Note that genotypic won't have an SD_value because it's not a continuous data type
        response_vals = panel_entry['BCG response characteristic']
        response_vals_list = response_vals.split('/')  # for when there's more than one possible response value
        HRSD_resp_list = HRSD_resp.split('/')  # type: ignore
        for response_val, HRSD in zip(response_vals_list, HRSD_resp_list):
            # log.write(response_val)
            # log.write(HRSD)
            response_allele1, response_allele2 = response_val.split('|')
            if result_allele1 == response_allele1 and result_allele2 == response_allele2:
                score = np.float64(HRSD)
                scored = True
            elif result_allele1 == response_allele2 and result_allele2 == response_allele1:
                score = np.float64(HRSD)
                scored = True
        
            if scored:
                log.write(f'score found for {panel_id}')
                break
            
    elif scoring_type == 'binary':
        if panel_entry[BIOMARKER_TYPE] == 'mod':
            beta_result = float64(result_entry['Result'])
            if beta_result < 0.8:  # TODO: check against classifier results
                binary_result = 'unmethylated'
                result = beta_result
            else:
                binary_result = 'methylated'
                result = beta_result
            control_val = panel_entry['Control Group or Baseline']
            if not np.isnan(control_val):
                raise ValueError("Haven't coded this yet")
            response_val = panel_entry['BCG response characteristic']
            nonresponse_val = panel_entry['BCG non response characteristic']
            if binary_result == response_val:
                score = float64(panel_entry['HR/SD Response'])
            elif binary_result == nonresponse_val:
                score = float64(panel_entry['HR/SD Non-Response'])
            else:
                raise ValueError("Shouldn't get here")
            scored = True
        else:
            raise ValueError("Haven't coded this yet")
    
    if scored:
        score_results_dict[panel_id] = (score, result, response, nonresponse)
        if response:
            log.write(f'Response score: {score}')
        if nonresponse:
            log.write(f'NonResponse score: {score}')
            
        # Get min/max scores
        min_score, max_score = get_min_max_scores(panel_entry)
        # Sense check
        if min_score > max_score:
            print(f'min score: {min_score}; max score: {max_score}')
            raise ValueError('Min score is bigger than max score')
            
        # print(f'min score: {min_score}; max score: {max_score}')
        min_scores_list.append(min_score)
        max_scores_list.append(max_score)
        
    else:
        # Check if alleles are present in the possible variants listed
        in_options = False
        option_list = result_entry['Result Options'].split('/')
        for option_vals in option_list:
            option_allele1, option_allele2 = option_vals.split('|')
            if result_allele1 == option_allele1 and result_allele2 == option_allele2:
                in_options = True
            elif result_allele1 == option_allele2 and result_allele2 == option_allele1:
                in_options = True
        
        if in_options:
            raise ValueError('Where did the score go?')
        else:
            log.write('Result alleles are not one of the options.')

score_only_results_dict = dict()
# results_dict = dict()
for panel_id, (score, result, response, nonresponse) in score_results_dict.items():
    # # if both response and nonresponse, leave score unchanged
    # if not nonresponse:
    #     score_only_results_dict[panel_id] = score * (-1)
    # else:
    #     score_only_results_dict[panel_id] = score
    score_only_results_dict[panel_id] = [result, np.float64(score)]
    # results_dict[panel_id] = result
        

total_score = np.nansum([v[1] for v in score_only_results_dict.values()])

log.write(str(total_score))
print(total_score)  # is the final number!! For now...

min_total = sum(min_scores_list)
max_total = sum(max_scores_list)

score_only_results_dict['Total (raw)'] = ['', total_score]
normalised = ((total_score - min_total) / (max_total - min_total))*100
score_only_results_dict['Total (normalised to 0.0 - 100.0)'] = ['', normalised]

score_results_df = pd.DataFrame.from_dict(score_only_results_dict, orient='index', columns = ['Result', 'Score'])
# score_results_df.reset_index(names = 'Panel ID', inplace=True)

score_results_df.to_csv(scores_csv)
    
log.close()
