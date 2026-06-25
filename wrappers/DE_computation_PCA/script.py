#############################################################
# wrapper for rule: normalize_and_visualize
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: normalize_and_visualize \n##\n")
f.close()

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/normalize_and_visualize.R " + \
          snakemake.input.count_data_original + " " + \
          snakemake.input.txi + " " + \
          snakemake.input.experiment_design + " " + \
          snakemake.wildcards.analysis_type + " " + \
          snakemake.params.condition_to_compare_vec + " " + \
          snakemake.output.output_dir + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)
