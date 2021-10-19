#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=28G,h_vmem=30G,h_fsize=200G
#$ -N step7-Rcounts-zandiHyde_Bipolar.LIBD
#$ -o ./logs/Rcounts-zandiHyde_Bipolar.txt
#$ -e ./logs/Rcounts-zandiHyde_Bipolar.txt
#$ -hold_jid pipeline_setup,step4-featCounts-zandiHyde_Bipolar.LIBD,step6-txQuant-zandiHyde_Bipolar.LIBD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/.create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data -e zandiHyde_Bipolar -p LIBD -l TRUE -c TRUE -t 5 -s reverse

echo "**** Job ends ****"
date
