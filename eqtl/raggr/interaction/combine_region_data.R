####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")

### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"

## load SNP data
load("../../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
snp = snp[-snpInd,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 


######################
# statistical model ##
######################
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))

mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]


######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")


#######################
####### do PCA ########
#######################

geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

rpkmA = rbind(geneRpkm,exonRpkm,jxnRp10m,txTpm)
snpA = snp
rownames(snpA) = snpMap$pos_hg19
mdsA = mds
pdA = pd[,c("BrainRegion","SAMPLE_ID","BrNum","RNum","PrimaryDx",
			"AgeDeath","Sex","Race","overallMapRate","totalAssignedGene")]

# pcaExp = prcomp(t(log2(rpkmA+1)))
# # kExp = num.sv(log2(rpkm+1), mod)
# # ExpPCs = pcaExp$x[,1:kExp]

# save(pcaExp, file="amyg_expressionPCs.rda")

geneMapA = as.data.frame(rowRanges(rse_gene))
exonMapA = as.data.frame(rowRanges(rse_exon))
jxnMapA = as.data.frame(rowRanges(rse_jxn))
txMapA = as.data.frame(rowRanges(rse_tx))
rownames(exonMapA) = exonMapA$exonPos


##################################################################################
##################################################################################
## sACC


## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")

### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"

## load SNP data
load("../../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
snp = snp[-snpInd,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 

######################
# statistical model ##
######################
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))

mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]

######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")

#######################
####### do PCA ########
#######################

geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

rpkmS = rbind(geneRpkm,exonRpkm,jxnRp10m,txTpm)
snpS = snp
rownames(snpS) = snpMap$pos_hg19
mdsS = mds
pdS = pd[,c("BrainRegion","SAMPLE_ID","BrNum","RNum","PrimaryDx",
			"AgeDeath","Sex","Race","overallMapRate","totalAssignedGene")]

# pcaExp = prcomp(t(log2(rpkmS+1)))
# # kExp = num.sv(log2(rpkm+1), mod)
# # ExpPCs = pcaExp$x[,1:kExp]

# save(pcaExp, file="sacc_expressionPCs.rda")




##################################################################################
##################################################################################
## DLPFC


## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pd = colData(rse_gene)

## load SNP data ## same 10777 snps used in other regions
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../rAggr_results_881.csv", stringsAsFactors=FALSE)
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 


######################
# statistical model ##
######################

pd$Dx = factor(pd$Dx,
	levels = c("Control", "Bipolar"))
	
mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]

######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")

#######################
####### do PCA ########
#######################

geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

rpkmD = rbind(geneRpkm,exonRpkm,jxnRp10m,txTpm)
snpD = snp
rownames(snpD) = snpMap$pos_hg19
mdsD = mds
pdD = pd[,c("Region","SAMPLE_ID","BrNum","RNum","Dx",
			"Age","Sex","Race","overallMapRate","totalAssignedGene")]
## match other regions	
pdD$SAMPLE_ID = unlist(lapply(pdD$SAMPLE_ID, function(x) paste(x, collapse=",")))
pdD$overallMapRate = unlist(lapply(pdD$overallMapRate, mean))	
pdD$totalAssignedGene = unlist(lapply(pdD$totalAssignedGene, mean))		
colnames(pdD) = c("BrainRegion","SAMPLE_ID","BrNum","RNum","PrimaryDx",
			"AgeDeath","Sex","Race","overallMapRate","totalAssignedGene")


# pcaExp = prcomp(t(log2(rpkmD+1)))
# # kExp = num.sv(log2(rpkm+1), mod)
# # ExpPCs = pcaExp$x[,1:kExp]

# save(pcaExp, file="dlpfc_expressionPCs.rda")

geneMapD = as.data.frame(rowRanges(rse_gene))
exonMapD = as.data.frame(rowRanges(rse_exon))
jxnMapD = as.data.frame(rowRanges(rse_jxn))
txMapD = as.data.frame(rowRanges(rse_tx))
rownames(exonMapD) = exonMapD$exonPos


##################################################################################
##################################################################################
## combined

## other data
snpCombined = cbind(snpA,snpS,snpD)
mdsCombined = rbind(mdsA,mdsS,mdsD)
pdCombined = rbind(pdA, pdS, pdD)

## combine expression
## amyg and sacc are the same
rpkmCombined = cbind(rpkmA, rpkmS)

identical(rownames(geneMapA), rownames(geneMapD))
identical(rownames(exonMapA), rownames(exonMapD))
identical(rownames(txMapA), rownames(txMapD))

## exon names
rownames(rpkmCombined)[(nrow(geneMapA)+1):(nrow(geneMapA)+nrow(exonMapA))] = rownames(exonMapA)
rownames(rpkmD)[(nrow(geneMapD)+1):(nrow(geneMapD)+nrow(exonMapD))] = rownames(exonMapD)

## jxn overlaps
rownames(rpkmCombined) = gsub("\\(\\-","\\(\\*", rownames(rpkmCombined))
rownames(rpkmCombined) = gsub("\\(\\+","\\(\\*", rownames(rpkmCombined))

jxnMapA$jxnPos = paste0(jxnMapA$seqnames,":",jxnMapA$start,"-",jxnMapA$end,"(*)")
jxnMapD$jxnPos = paste0(jxnMapD$seqnames,":",jxnMapD$start,"-",jxnMapD$end,"(*)")
length(which(!jxnMapA$jxnPos %in% jxnMapD$jxnPos))
length(which(!jxnMapD$jxnPos %in% jxnMapA$jxnPos))  ## nothing to drop

## put amyg and sacc in order of dlpfc
rpkmCombined = rpkmCombined[rownames(rpkmD),]

## combine
rpkmCombined = cbind(rpkmCombined, rpkmD)


save(snpCombined,mdsCombined,pdCombined,rpkmCombined,
		geneMapD, exonMapD, jxnMapD, txMapD, file="3regions_combined_data.rda")




### PCA on all 3 regions
pcaExp = prcomp(t(log2(rpkmCombined+1)))

save(pcaExp, file="3regions_expressionPCs.rda")







