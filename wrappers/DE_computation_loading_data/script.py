#############################################################
# wrapper for rule: load_count_data
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: load_count_data \n##\n")
f.close()

# creation of the TAB delimited table directly from the panda sample table
snakemake.params.sample_tab.to_csv(snakemake.output.experiment_design, index=False,sep="\t")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/load_count_data.R " + \
          snakemake.output.experiment_design + " " +\
          snakemake.input.count_tab + " " +\
          ",".join(snakemake.params.comparison_dir_list) + " " +\
          snakemake.input.gtf + " " + \
          snakemake.wildcards.analysis_type + " " + \
          str(snakemake.params.paired_replicates) + " " +\
          str(snakemake.params.use_custom_batch_effect_grouping) + " " +\
          str(snakemake.params.geneList) + " " + \
          str(snakemake.params.keepGene) + " " + \
          str(snakemake.params.chrmList) + " " + \
          str(snakemake.params.keepChrm) + " " + \
          str(snakemake.params.remove_genes_with_sum_read_count_threshold) + " " + \
          str(snakemake.params.remove_genes_with_mean_read_count_threshold) + " " + \
          snakemake.output.output_dir + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
shell(command)
