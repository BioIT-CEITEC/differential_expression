#############################################################
# wrapper for rule: DE_computation
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: DE_report \n##\n")
f.close()


if os.path.exists("DE_RSEM/all_condition_results/sample_to_sample_PCA_3D_batch.html"):
     command = "cp DE_RSEM/all_condition_results/sample_to_sample_PCA_3D_batch.html " + snakemake.output.html
elif os.path.exists("DE_RSEM/all_condition_results/sample_to_sample_PCA_3D.html"):
    command = "cp DE_RSEM/all_condition_results/sample_to_sample_PCA_3D.html " + snakemake.output.html
elif os.path.exists("DE_feature_count/all_condition_results/sample_to_sample_PCA_3D.html"):
    command = "cp DE_feature_count/all_condition_results/sample_to_sample_PCA_3D.html " + snakemake.output.html
elif os.path.exists("DE_feature_count/all_condition_results/sample_to_sample_PCA_3D.html"):
    command = "cp DE_feature_count/all_condition_results/sample_to_sample_PCA_3D.html " + snakemake.output.html
else:
    command = "touch " + snakemake.output.html

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)


command = " cp '"+os.path.abspath(os.path.dirname(__file__))+"/DE_report.Rmd' DE_report/"

shell(command)

command = """ Rscript -e "rmarkdown::render('DE_report/DE_report.Rmd', params=list(config = '""" + snakemake.params.config + """'))" """ +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " ls " + snakemake.output.html

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " rm DE_report/DE_report.Rmd"

shell(command)

command = " mv DE_report/DE_report.html " + snakemake.output.html

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
