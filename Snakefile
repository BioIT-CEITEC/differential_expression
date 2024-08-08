import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

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

#set analysis selected analysis types from config and rise exception if no selected
analysis = []
if config["featureCount"]:
    count_over_list = config['count_over'].split(",")
    if ("exon" in count_over_list):
        config["featureCount_exon"] = True
        analysis.append("featureCount_exon")
    if ("gene" in count_over_list):
        config["featureCount_gene"] = True
        analysis.append("featureCount_gene")
    if ("transcript" in count_over_list):
        config["featureCount_transcript"] = True
        analysis.append("featureCount_transcript")
    if ("three_prime_UTR" in count_over_list):
        config["featureCount_3pUTR"] = True
        analysis.append("featureCount_3pUTR")
    if ("five_prime_UTR" in count_over_list):
        config["featureCount_5pUTR"] = True
        analysis.append("featureCount_5pUTR")

if config["RSEM"]:
    analysis.append("RSEM")
if config["salmon_align"]:
    analysis.append("salmon_align")
if config["salmon_map"]:
    analysis.append("salmon_map")
if config["kallisto"]:
    analysis.append("kallisto")
if config["is_mirna"]:
    analysis.append("mirbase_canonical")
if len(analysis) == 0:
    raise ValueError("There was no RSEM or featureCount used in previous analysis!")


biotype_dir_list = config['biotypes'].split(",")

config["analysis_type"] = "|".join(analysis)
config["biotype_list"] = "|".join(biotype_dir_list)
config["comparison"] = "|".join(comparison_dir_list)

wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     analysis_type= "featureCount_exon|featureCount_gene|featureCount_transcript|featureCount_3pUTRn|featureCount_5pUTR|RSEM|salmon_map|salmon_align|kallisto|mirbase_canonical",
     #data_type= "tsv|RData"

os.makedirs("DE_report",exist_ok=True)

f=open("DE_report/DE_report.json", "w")
json.dump(config, f, indent=4)
f.close()

##### Target rules #####

def input_all(wildcards):
    input = {}
    if config["featureCount_exon"]:
        input["featureCount_exon"] = "DE_featureCount_exon/final_report.html"
    if config["featureCount_gene"]:
        input["featureCount_gene"] = "DE_featureCount_gene/final_report.html"
    if config["featureCount_transcript"]:
        input["featureCount_transcript"] = "DE_featureCount_transcript/final_report.html"
    if config["featureCount_3pUTR"]:
        input["featureCount_3pUTR"] = "DE_featureCount_3pUTR/final_report.html"
    if config["featureCount_5pUTR"]:
        input["featureCount_5pUTR"] = "DE_featureCount_5pUTR/final_report.html"
    if config["RSEM"]:
        input["RSEM"] = "DE_RSEM/final_report.html"
    if config["salmon_map"]:
        input["salmon_map"] = "DE_salmon_map/final_report.html"
    if config["salmon_align"]:
        input["salmon_align"] = "DE_salmon_align/final_report.html"
    if config["kallisto"]:
        input["kallisto"] = "DE_kallisto/final_report.html"
    if config["is_mirna"]:
        input["mirbase"] = "DE_mirbase_canonical/final_report.html"
    return input

rule all:
    input: expand("final_report.html")


##### Modules #####

include: "rules/differential_expression.smk"
