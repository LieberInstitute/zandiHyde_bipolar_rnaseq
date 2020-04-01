#!/bin/bash
#$ -cwd
#$ -l mem_free=80G,h_vmem=80G,h_fsize=200G
#$ -N plots
#$ -o logs/plots.txt
#$ -e logs/plots.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/conditional/plot_transcripts.R

echo "**** Job ends ****"
date
