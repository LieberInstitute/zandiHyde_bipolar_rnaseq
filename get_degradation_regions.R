###

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

## read manifest
load("data/zandiHypde_bipolar_rseGene_n511.rda")

## split by region
rse_gene_Amygdala = rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sACC = rse_gene[,rse_gene$BrainRegion == "sACC"]

## import degradation regions
bed_amyg = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Amygdala_RiboZero/bed/Amygdala_RibZero_degradation_top1000.bed")
bed_sacc = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/sACC_RiboZero/bed/sACC_RibZero_degradation_regions_top1000.bed")
bed_overall = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/bed/sACC_Plus_Amygdala_RibZero_degradation_top1000.bed")

#######################
##### Amygdala ###########
#######################

## designate bigwigs
forwardBw_amyg = paste0("preprocessed_data/Coverage/",
	rse_gene_Amygdala$SAMPLE_ID,".Forward.bw")
reverseBw_amyg = paste0("preprocessed_data/Coverage/",
	rse_gene_Amygdala$SAMPLE_ID, ".Reverse.bw")
all(file.exists(c(forwardBw_amyg,reverseBw_amyg))) # TRUE
names(forwardBw_amyg) = names(reverseBw_amyg) = rse_gene_Amygdala$SAMPLE_ID

## try coverage tool
covForward_amyg = coverage_bwtool(forwardBw_amyg, bed_amyg, strand = "+", 
	sumsdir = "degradation_amyg", bpparam = MulticoreParam(8))
covForward_amyg$bigwig_path = NULL
covForward_amyg$bigwig_file = NULL

covReverse_amyg = coverage_bwtool(reverseBw_amyg, bed_amyg, strand = "-", 
	sumsdir = "degradation_amyg", bpparam = MulticoreParam(8))
covReverse_amyg$bigwig_path = NULL
covReverse_amyg$bigwig_file = NULL

## combine
cov_rse_amyg = rbind(covForward_amyg, covReverse_amyg)	
rownames(cov_rse_amyg) = rowData(cov_rse_amyg)$name
cov_rse_amyg = cov_rse_amyg[bed_amyg$name,]

## divide by read length
assays(cov_rse_amyg)$counts = assays(cov_rse_amyg)$counts/100 # divide by read length

## make positive
assays(cov_rse_amyg)$counts = abs(assays(cov_rse_amyg)$counts) 

## save to final people
colData(cov_rse_amyg) = colData(rse_gene_Amygdala)
save(cov_rse_amyg, file = "data/degradation_rse_BipSeq_Amygdala.rda")

#######################
##### sACC ###########
#######################

## designate bigwigs
forwardBw_sacc = paste0("preprocessed_data/Coverage/",
	rse_gene_sACC$SAMPLE_ID,".Forward.bw")
reverseBw_sacc = paste0("preprocessed_data/Coverage/",
	rse_gene_sACC$SAMPLE_ID, ".Reverse.bw")
all(file.exists(c(forwardBw_sacc,reverseBw_sacc))) # TRUE
names(forwardBw_sacc) = names(reverseBw_sacc) = rse_gene_sACC$SAMPLE_ID

## try coverage tool
covForward_sacc = coverage_bwtool(forwardBw_sacc, bed_sacc, strand = "+", 
	sumsdir = "degradation_sacc", bpparam = MulticoreParam(8))
covForward_sacc$bigwig_path = NULL
covForward_sacc$bigwig_file = NULL

covReverse_sacc = coverage_bwtool(reverseBw_sacc, bed_sacc, strand = "-", 
	sumsdir = "degradation_sacc", bpparam = MulticoreParam(8))
covReverse_sacc$bigwig_path = NULL
covReverse_sacc$bigwig_file = NULL

## combine
cov_rse_sacc = rbind(covForward_sacc, covReverse_sacc)	
rownames(cov_rse_sacc) = rowData(cov_rse_sacc)$name
cov_rse_sacc = cov_rse_sacc[bed_sacc$name,]

## divide by read length
assays(cov_rse_sacc)$counts = assays(cov_rse_sacc)$counts/100 # divide by read length

## make positive
assays(cov_rse_sacc)$counts = abs(assays(cov_rse_sacc)$counts) 

## save to final people
colData(cov_rse_sacc) = colData(rse_gene_sACC)
save(cov_rse_sacc, file = "data/degradation_rse_BipSeq_sACC.rda")

#######################
##### Combined ###########
#######################

## designate bigwigs
forwardBw = paste0("preprocessed_data/Coverage/",
	rse_gene$SAMPLE_ID,".Forward.bw")
reverseBw = paste0("preprocessed_data/Coverage/",
	rse_gene$SAMPLE_ID, ".Reverse.bw")
all(file.exists(c(forwardBw,reverseBw))) # TRUE
names(forwardBw) = names(reverseBw) = rse_gene$SAMPLE_ID

## try coverage tool
covForward = coverage_bwtool(forwardBw, bed_overall, strand = "+", 
	sumsdir = "degradation_joint", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bed_overall, strand = "-", 
	sumsdir = "degradation_joint", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bed_overall$name,]

## divide by read length
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## save to final people
colData(cov_rse) = colData(rse_gene)
save(cov_rse, file = "data/degradation_rse_BipSeq_BothRegions.rda")
