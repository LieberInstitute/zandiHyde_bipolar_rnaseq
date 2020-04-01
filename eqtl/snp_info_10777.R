##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)

################
## load SNP data
load("../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## risk loci from PGC paper
indexLoci = read.csv("../PGC_risk_loci.csv", stringsAsFactors=FALSE) ## 881
indexIndex = which(snpMap$pos_hg19 %in% indexLoci$hg19POS)	# keep 456

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 10,777 snps from snpMap
snpMap = snpMap[keepIndex,]
keepIndex = which(riskLoci$hg19POS2 %in% snpMap$pos_hg19)	# keep 10,754 snps from riskLoci
riskLoci = riskLoci[keepIndex,]

snpMap$Status = ifelse(snpMap$pos_hg19 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status1 = ifelse(riskLoci$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status2 = ifelse(riskLoci$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")

## Also keep track of full list before dropping
riskLoci_full$hg19POS1 = paste0(riskLoci_full$SNP1_Chr, ":", riskLoci_full$SNP1_Pos) 
riskLoci_full$hg19POS2 = paste0(riskLoci_full$SNP2_Chr, ":", riskLoci_full$SNP2_Pos) 
riskLoci_full$SNP2_missing = "missing"
riskLoci_full$SNP2_missing[keepIndex] = "analyzed"
riskLoci_full$Status1 = ifelse(riskLoci_full$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci_full$Status2 = ifelse(riskLoci_full$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")



## the associated index SNP
## the associated index SNPs position in hg38 and h19
## and the distance and R2 of the SNP with the index SNP?
##  an indicator if the snp is in one of the PGC significant (n=31) loci

proxInd = match(snpMap$pos_hg19, riskLoci$hg19POS2)
snpMap$IndexSNP = riskLoci$SNP1_Name[proxInd]
snpMap$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]
snpMap$IndexSNP_hg38POS = snpMap$pos_hg38[match(snpMap$IndexSNP_hg19POS,snpMap$pos_hg19)]
snpMap$Distance = riskLoci$Distance[proxInd]
snpMap$R_squared = riskLoci$R_squared[proxInd]

## keep only genomewide
indexLoci = read.csv("../PGC_risk_loci.csv", stringsAsFactors=FALSE) ## 881
indexLoci = indexLoci[indexLoci$genomewide,]

snpMap$IndexSNP_genomewide = snpMap$IndexSNP_hg19POS %in% indexLoci$hg19POS




write.csv(snpMap, file="snpMap_10777.csv")



