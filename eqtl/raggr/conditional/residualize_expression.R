##
library(jaffelab)
library(SummarizedExperiment)

################### amygdala

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
pd = colData(rse_gene)

## load SNP data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[pd$BrNum,]
rownames(mds) = pd$RNum

load("../rdas/pcs_4features_amyg.rda", verbose=TRUE)
##
modG = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + genePCs, data = pd)
modE = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + exonPCs, data = pd)
modJ = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + jxnPCs, data = pd)
modT = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + txPCs, data = pd)

################
## load expression
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, modG, P=1)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, modE, P=1)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, modJ, P=1)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, modT, P=1)

amyg_log2exprs_clean = rbind(gExprs,eExprs,jExprs,tExprs)
pd_amyg = pd



################### sacc

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
pd = colData(rse_gene)

## load SNP data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[pd$BrNum,]
rownames(mds) = pd$RNum

load("../rdas/pcs_4features_sacc.rda", verbose=TRUE)
##
modG = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + genePCs, data = pd)
modE = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + exonPCs, data = pd)
modJ = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + jxnPCs, data = pd)
modT = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + txPCs, data = pd)

################
## load expression
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, modG, P=1)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, modE, P=1)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, modJ, P=1)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, modT, P=1)

sacc_log2exprs_clean = rbind(gExprs,eExprs,jExprs,tExprs)
pd_sacc = pd


################### dlpfc

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pd = colData(rse_gene)

## load SNP data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
mds = mds[pd$BrNum,]
rownames(mds) = pd$RNum

load("../rdas/pcs_4features_dlpfc.rda", verbose=TRUE)
##
modG = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + genePCs, data = pd)
modE = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + exonPCs, data = pd)
modJ = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + jxnPCs, data = pd)
modT = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + txPCs, data = pd)

################
## load expression
## dlpfc
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, modG, P=1)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, modE, P=1)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, modJ, P=1)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, modT, P=1)

dlpfc_log2exprs_clean = rbind(gExprs,eExprs,jExprs,tExprs)
pd_dlpfc = pd


save(pd_amyg, amyg_log2exprs_clean, 
	 pd_sacc, sacc_log2exprs_clean, 
	 pd_dlpfc, dlpfc_log2exprs_clean, file="residualized_exprs_3regions.rda")







