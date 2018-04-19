#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=120G,h_vmem=120G,h_fsize=200G
#$ -N dlpfc_eQTL_bipolar
#$ -o logs/eQTL_dlpfc.txt
#$ -e logs/eQTL_dlpfc.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/raggr_run_eqtls_dlpfc.R

echo "**** Job ends ****"
date
