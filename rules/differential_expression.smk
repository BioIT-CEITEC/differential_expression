
import pandas as pd
from helper import define_variable





###########################################
# DEFINITION OF VARIABLES
#
cfg = pd.DataFrame(config)

REF_DIR = define_variable(cfg, "REF_DIR")
PROJECT_NAME = define_variable(cfg, "PROJECT_NAME")
PROJECT_DIR = define_variable(cfg, "PROJECT_DIR")
INPUTS_DIR = define_variable(cfg, "INPUTS_DIR")
ADIR = define_variable(cfg, "ANALYSIS_DIR")

if not "ref_from_trans_assembly" in cfg:
    cfg["ref_from_trans_assembly"] = "F"

# REF_DIR = "/mnt/ssd/ssd_3/references"
# PROJECT_NAME = cfg['project'].tolist()[0].replace("/",".")
# PROJECT_DIR = os.path.join(cfg['project_owner'].tolist()[0]
#                     ,"sequencing_results"
#                     ,"projects"
#                     ,cfg['project'].tolist()[0])
# if "test_new_staging" in cfg['project'].tolist()[0]:
#     PROJECT_DIR = "stage"+cfg['stage_id'].tolist()[0]+"_"+cfg['project'].tolist()[0].replace("/",".")
#
# INPUTS_DIR = os.path.join(PROJECT_DIR,"input_files")
# ADIR = os.path.join(PROJECT_DIR,cfg['analysis_name'].tolist()[0])


####################################
# FINAL RESULTING FILES
#


def final_variant_calling_report_input(wildcards):

    if cfg['conditions_to_compare'].tolist()[0] == "all":
        condition_list = sorted(list(set(cfg['condition'].tolist())))
    else:
        condition_list = cfg['conditions_to_compare'].tolist()[0].split(",")

    comparison_dir_list = list()
    for condition1 in condition_list:
        if ':' in condition1:
            conditions = condition1.split(":")
            comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
        else:
            for condition2 in condition_list[condition_list.index(condition1):]:
                if ':' not in condition2 and condition2 != condition1:
                    comparison_dir_list.append(condition2 + "_vs_" + condition1)

    biotype_dir_list = cfg['biotypes'].tolist()[0].split(",")
    return expand(ADIR + "/{comparison}/{biotype}/edgeR.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)




def DE_computation_input(wildcards):
    if cfg['analysis_subclass'].tolist()[0] == "featureCounts":
        return ADIR+"/complete.featureCounts.tsv"
    else:
        return ADIR+"/complete.RSEM.RData"


rule DE_computation:
    input:  expression_tab=DE_computation_input,
            cfg_tab = ADIR + "/" + PROJECT_NAME + "." + cfg['analysis_class'].tolist()[0] + ".config.json",
            biotype_groups = REF_DIR + "/general/default/annot/biotypes_list_mod.txt",
            sqlite = expand("{dir}/{organism}/{ref}/annot/{ref}.sqlite.gz", dir=REF_DIR, organism=cfg["organism"].tolist()[0], ref=cfg["reference"].tolist()[0]),
    output: table = ADIR + "/{comparison}/{biotype}/edgeR.tsv",
    params: count_type = cfg['analysis_subclass'].tolist()[0],
            organism = cfg["organism"].tolist()[0],
            use_tag_to_pair_samples = cfg["use_tag_to_pair_samples"].min(),
            ref_from_trans_assembly = cfg["ref_from_trans_assembly"].min()
    log:    run = ADIR + "/{comparison}/{biotype}/{comparison}.{biotype}.DE.run.log"
    conda:  "../wrappers/DE_computation/env.yaml"
    script: "../wrappers/DE_computation/script.py"

rule analysis_RSEM_table:
    input:  RSEM = set(expand(INPUTS_DIR+"/rsem_counts/{full_name}.genes.results",full_name = cfg['full_name'].tolist()))
    output: RSEM = ADIR+"/complete.RSEM.RData"
    params: gtf = expand("{dir}/{organism}/{ref}/annot/{ref}.gtf", dir=REF_DIR, organism=cfg["organism"].min(), ref=cfg["reference"].min()),
            ref_from_trans_assembly = cfg["ref_from_trans_assembly"].min()
    log:    run = ADIR + "/complete.RSEM.log"
    conda:  "../wrappers/analysis_RSEM_table/env.yaml"
    script: "../wrappers/analysis_RSEM_table/script.py"

rule analysis_featureCounts_table:
    input:  featureCounts = set(expand(INPUTS_DIR+"/feature_counts/{full_name}.featureCounts.tsv",full_name = cfg['full_name'].tolist()))
    output: table = ADIR+"/complete.featureCounts.tsv",
    log:    run = ADIR + "/complete.featureCounts.log"
    conda:  "../wrappers/project_featureCounts_table/env.yaml"
    script: "../wrappers/project_featureCounts_table/script.py"
