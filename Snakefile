import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

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


# if config["is_paired"] == False:
#     read_pair_tags = [""]
#     paired = "SE"
# else:
#     read_pair_tags = ["_R1","_R2"]
#     paired = "PE"

#set analysis selected analysis types from config and rise exception if no selected
selected_analysis_types = []
if config["feature_count"]:
    selected_analysis_types.append("feature_count")
if config["RSEM"]:
    selected_analysis_types.append("RSEM")

# if len(selected_analysis_types) == 0:
#     !!exception


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
