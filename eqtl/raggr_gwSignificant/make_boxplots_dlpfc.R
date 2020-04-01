##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pd = colData(rse_gene)

## load SNP data ## same 10777 snps used in other regions
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

# ## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
# snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
# snpMap = snpMap[-snpInd,]
# snp = snp[-snpInd,]

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 10777 snps from snpMap
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)

mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum



################
## load table
dlp = read.csv("raggr_31_snps_dlpfc_eqtls_fdr01.csv", row.names=1)
dlp$Symbol = as.character(dlp$Symbol)
dlp$SNP = as.character(dlp$SNP)
dlp$gene = as.character(dlp$gene)
dlp$IndexSNP = as.character(dlp$IndexSNP)

load("../raggr/rdas/pcs_4features_dlpfc.rda", verbose=TRUE)
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

exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)
exprsAdj2 = rbind(log2(geneRpkm+1),log2(exonRpkm+1),log2(jxnRp10m+1),log2(txTpm+1))

dlp2 = dlp[order(dlp$pvalue, decreasing=FALSE),]

pdf("dlpfc_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:1000) {
	symi = dlp2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlp2[i,"SNP"]
	feati = dlp2[i,"gene"]
	p_i = signif(dlp2[i,"pvalue"],3)
	typei = dlp2[i,"Type.1"]

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



#### plot proxy and index next to each other
pdf("dlpfc_top_eqtl_adj_index.pdf", h=8, w=8)
par(mfrow=c(2,2), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:200) {
	symi = dlp2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlp2[i,"SNP"]
	snpicoord = dlp2[i,"hg19POS"]
	indexi = dlp2[i,"IndexSNP"]
	indexicoord = dlp2[i,"IndexSNP_hg19POS"]
	feati = dlp2[i,"gene"]
	pindex_i = signif(dlp2$pvalue[which(dlp2$SNP==indexi & dlp2$gene==feati)],3)
	p_i = signif(dlp2[i,"pvalue"],3)
	typei = dlp2[i,"Type.1"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=paste0(snpi,"\n",snpicoord), ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
	
	boxplot(exprsAdj[feati,] ~ snp[indexi,],
			xlab=paste0(indexi,"\n",indexicoord), ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[indexi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[indexi,])+2,cex=1.5)			   
	legend("top", paste0("p=",pindex_i))
}
dev.off()



