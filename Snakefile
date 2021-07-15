import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

#GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_REF_PATH = "/mnt/ssd/ssd_3/references"

# DNA parameteres processing
#

if not "ref_from_trans_assembly" in config:
    config["ref_from_trans_assembly"] = "F"




# # Reference processing
# #
# if config["lib_ROI"] != "wgs":
#     # setting reference from lib_ROI
#     f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","DNA_ROI.json"))
#     lib_ROI_dict = json.load(f)
#     f.close()
#     config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
#
#
# # setting organism from reference
# f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","genomic_references.json"),)
# reference_dict = json.load(f)
# f.close()
# config["organism"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]


##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+"

##### Target rules #####

rule all:
    input:  expand("results/{comparison}/{biotype}/edgeR.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)


##### Modules #####

include: "rules/differential_expression.smk"
