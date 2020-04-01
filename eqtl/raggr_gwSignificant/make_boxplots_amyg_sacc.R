##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")

### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"

## load SNP data
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
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
riskLoci = read.csv("../raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 



################
## load table
amyg = read.csv("raggr_31_snps_amyg_eqtls_fdr01.csv", row.names=1)
sacc = read.csv("raggr_31_snps_sacc_eqtls_fdr01.csv", row.names=1)
dlp = read.csv("raggr_31_snps_dlpfc_eqtls_fdr01.csv", row.names=1)


##
pd$Dx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))
mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)


################
## load expression
## amyg
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, mod, P=2)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, mod, P=2)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, mod, P=2)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, mod, P=2)


exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)
amyg$Symbol = as.character(amyg$Symbol)

pdf("amyg_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:100) {
	symi = amyg[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = amyg[i,"SNP"]
	feati = amyg[i,"gene"]
	p_i = signif(amyg[i,"pvalue"],3)
	typei = amyg[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
dev.off()










