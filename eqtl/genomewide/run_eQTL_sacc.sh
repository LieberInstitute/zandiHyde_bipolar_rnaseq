#!/bin/bash
#$ -cwd
#$ -l mem_free=150G,h_vmem=150G,h_fsize=200G
#$ -N sacc_eQTL_bipolar
#$ -o logs/eQTL_sacc.txt
#$ -e logs/eQTL_sacc.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/genomewide/genomewide_run_eqtls_sacc.R

echo "**** Job ends ****"
date
