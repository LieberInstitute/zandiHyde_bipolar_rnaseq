#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=1G,h_vmem=1G
#$ -j y
#$ -o logs/clean_trash_sacc.txt
#$ -m e

echo "**** Job starts ****"
date

rm -r /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/sacc_gene/tmp_files
rm -r /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/sacc_gene/out_files

echo "**** Job ends ****"
date
