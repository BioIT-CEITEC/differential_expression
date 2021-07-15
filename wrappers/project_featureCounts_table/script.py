#############################################################
# wrapper for rule: analysis_featureCounts_table
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: analysis_featureCounts_table \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_featureCount_tabs.R "+\
            snakemake.output.table+ " " +\
            " ".join(snakemake.input.featureCounts)

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
