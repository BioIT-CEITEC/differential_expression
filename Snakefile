import os
import pandas as pd
import json
from snakemake.utils import min_version
import re

min_version("5.18.0")

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")

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
    condition_list = sorted(sample_tab.condition.unique())
    condition_list_first = [condition for condition in condition_list if
                            not re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list_second = [condition for condition in condition_list if
                             re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list = condition_list_first + condition_list_second
    comparison_dir_list = get_comparison_dir_list(condition_list)
else:
    comparison_dir_list = get_comparison_dir_list(config['conditions_to_compare'].split(","))
    condition_list = set(config['conditions_to_compare'].replace(':',',').split(","))
    if not config['keep_not_compared_samples_for_normalization']:
        sample_tab = sample_tab[sample_tab['condition'].isin(condition_list)]

#
analysis = []
if config["feature_count"]:
    analysis.append("feature_count")
if config["RSEM"]:
    analysis.append("RSEM")
if config["salmon_align"]:
    analysis.append("salmon_align")
if config["salmon_map"]:
    analysis.append("salmon_map")
if config["kallisto"]:
    analysis.append("kallisto")

biotype_dir_list = config['biotypes'].split(",")

config["analysis_type"] = "|".join(analysis)
config["biotype_list"] = "|".join(biotype_dir_list)
config["comparison"] = "|".join(comparison_dir_list)

#set analysis selected analysis types from config and rise exception if no selected
selected_analysis_types = []
if config["feature_count"]:
    selected_analysis_types.append("feature_count")
if config["RSEM"]:
    selected_analysis_types.append("RSEM")
if config["salmon_map"]:
    selected_analysis_types.append("salmon_map")
if config["salmon_align"]:
    selected_analysis_types.append("salmon_align")
if config["kallisto"]:
    selected_analysis_types.append("kallisto")

if len(selected_analysis_types) == 0:
    raise ValueError("There was no RSEM or feature_count used in previous analysis!")


wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     analysis_type= "feature_count|RSEM|salmon_map|salmon_align|kallisto",
     #data_type= "tsv|RData"

os.makedirs("DE_report",exist_ok=True)

f=open("DE_report/DE_report.json", "w")
json.dump(config, f, indent=4)
f.close()

##### Target rules #####

def input_all(wildcards):
    input = {}
    if config["feature_count"]:
        input["feature_count"] = "DE_feature_count/final_report.html"
    if config["RSEM"]:
        input["RSEM"] = "DE_RSEM/final_report.html"
    if config["salmon_map"]:
        input["salmon_map"] = "DE_salmon_map/final_report.html"
    if config["salmon_align"]:
        input["salmon_align"] = "DE_salmon_align/final_report.html"
    if config["kallisto"]:
        input["kallisto"] = "DE_kallisto/final_report.html"
    return input

rule all:
    input: expand("final_report.html")


##### Modules #####

include: "rules/differential_expression.smk"
