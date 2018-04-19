#!/bin/bash
#$ -cwd
#$ -l mem_free=120G,h_vmem=120G,h_fsize=200G
#$ -N amyg_eQTL_bipolar
#$ -o logs/eQTL_amyg.txt
#$ -e logs/eQTL_amyg.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/genomewide/genomewide_run_eqtls_amyg.R

echo "**** Job ends ****"
date
