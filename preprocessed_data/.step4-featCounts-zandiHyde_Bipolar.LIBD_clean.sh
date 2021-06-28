#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-zandiHyde_Bipolar.LIBD_clean
#$ -o ./logs/featCounts-zandiHyde_Bipolar_clean.txt
#$ -e ./logs/featCounts-zandiHyde_Bipolar_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-zandiHyde_Bipolar.LIBD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/Counts/junction/tmpdir

echo "**** Job ends ****"
date
