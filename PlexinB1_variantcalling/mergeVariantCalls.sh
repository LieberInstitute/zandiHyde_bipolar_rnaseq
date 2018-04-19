#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=10G
#$ -N mergeVariants-Bipolar
#$ -o ./logs/mergeVariants-Bipolar.txt
#$ -e ./logs/mergeVariants-Bipolar.txt
#$ -step8-callVariants-Bipolar
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/PlexinB1_variantcalling/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/PlexinB1_variantcalling/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date

