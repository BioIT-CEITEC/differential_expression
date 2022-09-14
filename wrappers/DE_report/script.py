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
