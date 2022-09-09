import os
import pandas as pd
import json
from snakemake.utils import min_version
import re

min_version("5.18.0")

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if (sample_tab.condition == "").all():
    raise ValueError("There are no conditions set for samples!")


def get_comparison_dir_list(condition_list):
    comparison_dir_list = list()
    for condition1 in condition_list:
        if ':' in condition1:
            conditions = condition1.split(":")
            comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
        else:
            for condition2 in condition_list[condition_list.index(condition1):]:
                if ':' not in condition2 and condition2 != condition1:
                    comparison_dir_list.append(condition1 + "_vs_" + condition2)
    return comparison_dir_list


## create list of conditions
if config['conditions_to_compare'] == "all":
    condition_list = sorted(sample_tab.condition)
    condition_list_first = [condition for condition in condition_list if
                            re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list_second = [condition for condition in condition_list if
                             not re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list = condition_list_first + condition_list_second
    comparison_dir_list = get_comparison_dir_list(condition_list)
else:
    comparison_dir_list = get_comparison_dir_list(config['conditions_to_compare'].split(","))
    condition_list = set(config['conditions_to_compare'].replace(':',',').split(","))
    if not config['keep_not_compared_samples_for_normalization']:
        sample_tab = sample_tab[sample_tab['condition'].isin(condition_list)]



#set analysis selected analysis types from config and rise exception if no selected
selected_analysis_types = []
if config["feature_count"]:
    selected_analysis_types.append("feature_count")
if config["RSEM"]:
    selected_analysis_types.append("RSEM")

if len(selected_analysis_types) == 0:
    raise ValueError("There was no RSEM or feature_count used in previous analysis!")


wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     analysis_type= "feature_count|RSEM",
     #data_type= "tsv|RData"


##### Target rules #####

def input_all(wildcards):
    input = {}
    if config["feature_count"]:
        input["feature_count"] = "DE_feature_count/final_report.html"
    if config["RSEM"]:
        input["RSEM"] = "DE_RSEM/final_report.html"
    return input

rule all:
    input: expand("DE_{analysis_type}/final_report.html",analysis_type=selected_analysis_types)


##### Modules #####

include: "rules/differential_expression.smk"
