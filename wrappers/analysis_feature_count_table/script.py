#############################################################
# wrapper for rule: analysis_feature_count_table
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: analysis_feature_count_table \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_feature_count_tabs.R "+\
            snakemake.output.table+ " " +\
            " ".join(snakemake.input.feature_count)

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
