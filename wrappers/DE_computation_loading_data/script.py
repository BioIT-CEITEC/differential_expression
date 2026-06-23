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

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/load_count_data.R " + \
          snakemake.input.count_tab + " " + \
          snakemake.params.experiment_design + " " + \
          snakemake.input.gtf + " " + \
          snakemake.wildcards.analysis_type + " " + \
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
shell(command)
