

def final_variant_calling_report_input(wildcards):
    input = {}

    if (sample_tab.condition != "").any() or (sample_tab.tag != "").any():
        if config['conditions_to_compare'] == "all":
            condition_list = sorted(sample_tab.condition)
        else:
            condition_list = config['conditions_to_compare'].split(",")

        comparison_dir_list = list()
        for condition1 in condition_list:
            if ':' in condition1:
                conditions = condition1.split(":")
                comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
            else:
                for condition2 in condition_list[condition_list.index(condition1):]:
                    if ':' not in condition2 and condition2 != condition1:
                        comparison_dir_list.append(condition2 + "_vs_" + condition1)

        biotype_dir_list = config['biotypes'].split(",")

        input['tsv'] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/edgeR.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)

        if config["ref_from_trans_assembly"] != False:
            input['trans_ids_map'] = expand("{ref_dir}/annot/{ref}.transdecoder_ids_map", ref_dir=reference_directory,ref=config["reference"])[0]
    else:
        raise ValueError("There is no conditions or tag for samples!")
    return input


rule final_DE_report:
    input: unpack(final_variant_calling_report_input)
    output: html = "results/{analysis_type}_final_report.html"
    params: config = "config.json",
            paired = paired,
            count_type="{analysis_type}",
            contaminants="all",  ### přidat do configu!
    shell: "touch {output.html}"

def DE_computation_input(wildcards):
    input={}
    input["cfg_tab"] = "config.json"
    input["biotype_groups"] = os.path.join(GLOBAL_REF_PATH,"general/default/annot/biotypes_list_mod.txt")
    input["sqlite"] = expand("{ref_dir}/annot/{ref}.sqlite.gz",ref_dir=reference_directory,ref=config["reference"])[0]
    if wildcards.analysis_type == "feature_count":
        input["expression_tab"] = "results/analysis_{analysis_type}_table/complete.{analysis_type}.tsv"
    if wildcards.analysis_type == "RSEM":
        input["expression_tab"] = "results/analysis_{analysis_type}_table/complete.{analysis_type}.RData"
    return input


rule DE_computation:
    input: unpack(DE_computation_input)
    output: table = "results/DE_{analysis_type}/{comparison}/{biotype}/edgeR.tsv",
    params: count_type="{analysis_type}",
            organism = config["organism"],
            use_tag_to_pair_samples = config["use_tag_to_pair_samples"],
            ref_from_trans_assembly = config["ref_from_trans_assembly"]
    log:    "logs/DE_{analysis_type}/{comparison}/{biotype}/{comparison}.{biotype}.DE.log"
    conda:  "../wrappers/DE_computation/env.yaml"
    script: "../wrappers/DE_computation/script.py"


rule analysis_RSEM_table:
    input:  RSEM = expand("qc_reports/{sample}/RSEM/{sample}.genes.results",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: RSEM = "results/analysis_RSEM_table/complete.RSEM.RData"
    params: ref_from_trans_assembly = config["ref_from_trans_assembly"]
    log:    "logs/RSEM/complete.RSEM.log"
    conda:  "../wrappers/analysis_RSEM_table/env.yaml"
    script: "../wrappers/analysis_RSEM_table/script.py"

rule analysis_feature_count_table:
    input:  feature_count = expand("qc_reports/{sample}/feature_count/{sample}.feature_count.tsv",sample=sample_tab.sample_name)
    output: table = "results/analysis_feature_count_table/complete.feature_count.tsv",
    log:    "logs/feature_count/complete.feature_count.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"
