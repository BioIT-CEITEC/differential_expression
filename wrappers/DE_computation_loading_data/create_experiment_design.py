#############################################################
# wrapper for rule: create_experiment_design
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: create_experiment_design \n##\n")
f.close()

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/create_experiment_design.R " + \
          snakemake.input.sample_tab + " " + \
          snakemake.output.experiment_design + " " + \
          str(snakemake.params.paired_replicates).upper() + " " + \
          str(snakemake.params.use_custom_batch_effect_grouping).upper() + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)
