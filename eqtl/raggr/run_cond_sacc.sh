#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=200G
#$ -N sacc_conditional
#$ -o logs/conditional_sacc_snp.txt
#$ -e logs/conditional_sacc_snp.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/pgc_conditional_sacc.R

echo "**** Job ends ****"
date
