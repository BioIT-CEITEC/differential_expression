#############################################################
# wrapper for rule: analysis_RSEM_table
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: analysis_RSEM_table \n##\n")
f.close()

if not snakemake.params.ref_from_trans_assembly:
    gtf_file = "ref_not_from_trans_assembly"
else:
    gtf_file = snakemake.params.gtf[0]

command = "$CONDA_PREFIX/bin/Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_RSEM_tabs.R "+\
            snakemake.output.RSEM + " " +\
            gtf_file + " " +\
            " ".join(snakemake.input.RSEM) + " >> " + snakemake.log.run + " 2>&1 "

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
