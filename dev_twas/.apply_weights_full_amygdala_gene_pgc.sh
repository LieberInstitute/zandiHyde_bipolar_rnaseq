
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N apply_weights_full_amygdala_gene_pgc
#$ -o ./logs/apply_weights_full_amygdala_gene_pgc.txt
#$ -e ./logs/apply_weights_full_amygdala_gene_pgc.txt
#$ -m e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

## Load dependencies
module load fusion_twas/github
module load conda_R/4.0

## List current modules
module list

## Choose the correct GWAS summary statistics file
if [ "pgc" == "pgc" ]
then
    summstatsfile="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/PGC_BIP_hg38_clean.txt"
else
    echo "Unexpected pgc input"
fi

## Apply weights for the given region/feature pair and the given GWAS summary statistics
mkdir -p "amygdala_gene/pgc"


for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome ${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R     --sumstats ${summstatsfile}     --weights /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/amygdala_gene/amygdala_gene.pos     --weights_dir /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/amygdala_gene/     --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR.     --chr ${chr}     --out amygdala_gene/pgc/pgc.${chr}.dat

    echo ""
    echo "making plots for chromosome ${chr}"
    date
    echo ""

## companion post-processing step (plots only for genes)
if [ "gene" == "gene" ]
then
    Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R         --sumstats /${summstatsfile}         --input amygdala_gene/pgc/pgc.${chr}.dat         --out amygdala_gene/pgc/pgc.${chr}.analysis         --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR.         --chr ${chr}         --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr         --glist_path "/jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/glist-hg38"
else
    Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R         --sumstats /${summstatsfile}         --input amygdala_gene/pgc/pgc.${chr}.dat         --out amygdala_gene/pgc/pgc.${chr}.analysis         --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR.         --chr ${chr}         --locus_win 100000 --verbose 2 --plot_corr
fi

done

echo "**** Job ends ****"
date
