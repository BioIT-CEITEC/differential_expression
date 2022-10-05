import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

if not "ref_from_trans_assembly" in config:
    config["ref_from_trans_assembly"] = "F"


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


if config["is_paired"] == False:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"


wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     analysis_type= "feature_count|RSEM|salmon_map|salmon_aln|kallisto",
     #data_type= "tsv|RData"


##### Target rules #####

def input_all(wildcards):
    input = {}
    if config["feature_count"]:
        input["feature_count"] = "results/feature_count_final_report.html"
    if config["RSEM"]:
        input["RSEM"] = "results/RSEM_final_report.html"
    if config["salmon_map"]:
        input["salmon_map"] = "results/salmon_map_final_report.html"
    if config["salmon_align"]:
        input["salmon_aln"] = "results/salmon_aln_final_report.html"
    if config["kallisto"]:
        input["kallisto"] = "results/kallisto_final_report.html"
    return input

rule all:
    input:  unpack(input_all)


##### Modules #####

include: "rules/differential_expression.smk"
