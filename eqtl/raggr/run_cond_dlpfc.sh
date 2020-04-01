#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=200G
#$ -N dlpfc_conditional
#$ -o logs/conditional_dlpfc_snp.txt
#$ -e logs/conditional_dlpfc_snp.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/pgc_conditional_dlpfc.R

echo "**** Job ends ****"
date
