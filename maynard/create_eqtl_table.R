##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)

################
## load SNP data
load("../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

################
## load eQTLs

load("../eqtl/genomewide/mergedEqtl_output_sacc_genomewide_4features_FDR01.rda", verbose=TRUE)
sacc = allEqtlFDR01
load("../eqtl/genomewide/mergedEqtl_output_amyg_genomewide_4features_FDR01.rda", verbose=TRUE)
amyg = allEqtlFDR01

## add coords to eQTL snps
indexInd = match(amyg$snps, snpMap$SNP) ## row of snp in snpMap
amyg$snp_chr = snpMap$CHR[indexInd]
amyg$snp_pos_hg38 = snpMap$pos_hg38[indexInd]
amyg$snp_pos_hg19 = snpMap$pos_hg19[indexInd]
amyg = amyg[,c(1,12:14,2:11)]

indexInd = match(sacc$snps, snpMap$SNP)
sacc$snp_chr = snpMap$CHR[indexInd]
sacc$snp_pos_hg38 = snpMap$pos_hg38[indexInd]
sacc$snp_pos_hg19 = snpMap$pos_hg19[indexInd]
sacc = sacc[,c(1,12:14,2:11)]


################
## load SNP list

riskLoci = read.csv("1425669816_060618_rAggr_output_MDD_GWASsig_chrpos_queried.csv")
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 


################
## eQTLs from snp list

indSacc = which(sacc$snp_pos_hg19 %in% riskLoci$hg19POS2)
indAmyg = which(amyg$snp_pos_hg19 %in% riskLoci$hg19POS2)

sacc2 = sacc[indSacc,]
amyg2 = amyg[indAmyg,]

# summary table
snpHits = data.frame(snp=unique(c(sacc2$snps,amyg2$snps)))	
indexInd = match(snpHits$snp, snpMap$SNP) ## row of snp in snpMap
snpHits$snp_chr = snpMap$CHR[indexInd]
snpHits$snp_pos_hg38 = snpMap$pos_hg38[indexInd]
snpHits$snp_pos_hg19 = snpMap$pos_hg19[indexInd]
# T/F
snpHits$eQTLinSacc = snpHits$snp_pos_hg19 %in% sacc2$snp_pos_hg19
snpHits$eQTLinAmyg = snpHits$snp_pos_hg19 %in% amyg2$snp_pos_hg19

write.csv(snpHits,"maynard_raggr_hydebipolar_eQTLs.csv")


############################
####### add proxy info

## Sort riskLoci by R2 so highest linked are chosen (i.e. index matches with itself)
riskLoci = riskLoci[order(riskLoci$R_squared, decreasing=TRUE),]

## Is SNP index snp or proxy
sacc2$Status = ifelse(sacc2$snp_pos_hg19 %in% riskLoci$hg19POS1, "Index", "Proxy")
## What is the index snp for each row
proxInd = match(sacc2$snp_pos_hg19, riskLoci$hg19POS2)
sacc2$IndexSNP = riskLoci$SNP1_Name[proxInd]
sacc2$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]
# Match proxy+index row
proxInd = match(paste(sacc2$snp_pos_hg19,sacc2$IndexSNP_hg19POS), paste(riskLoci$hg19POS2,riskLoci$hg19POS1))
sacc2$r2IndexSNP_ProxySNP = riskLoci$R_squared[proxInd]

#############################

## Is SNP index snp or proxy
amyg2$Status = ifelse(amyg2$snp_pos_hg19 %in% riskLoci$hg19POS1, "Index", "Proxy")
## What is the index snp for each row
proxInd = match(amyg2$snp_pos_hg19, riskLoci$hg19POS2)
amyg2$IndexSNP = riskLoci$SNP1_Name[proxInd]
amyg2$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]
# Match proxy+index row
proxInd = match(paste(amyg2$snp_pos_hg19,amyg2$IndexSNP_hg19POS), paste(riskLoci$hg19POS2,riskLoci$hg19POS1))
amyg2$r2IndexSNP_ProxySNP = riskLoci$R_squared[proxInd]


sacc2 = sacc2[,c(1:4,15:18,5:14)]
amyg2 = amyg2[,c(1:4,15:18,5:14)]

sacc2$gencodeTx = paste(sacc2$gencodeTx, collapse = ',')
amyg2$gencodeTx = paste(amyg2$gencodeTx, collapse = ',')


write.csv(sacc2,"maynard_raggr_hydebipolar_eQTLs_sacc.csv")
write.csv(amyg2,"maynard_raggr_hydebipolar_eQTLs_amyg.csv")










