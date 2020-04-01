####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")

### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"

## load SNP data
load("rdas/overlappingSNPs.rda")  # snpMapKeep
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

snpInd = which(rownames(snpMap) %in% rownames(snpMapKeep) & !is.na(snpMap$pos_hg38))
snpMap = snpMap[snpInd,]
snp = snp[snpInd,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 


################
## load table
load("mergedEqtl_output_sacc_genomewide_4features_FDR01.rda", verbose=TRUE)
sacc = allEqtlFDR01


##
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))
mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)


################
## load expression
## sacc
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, mod, P=1)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, mod, P=1)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, mod, P=1)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, mod, P=1)


exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)
sacc$Symbol = as.character(sacc$Symbol)

saccG = sacc[which(sacc$Type=="Gene"),]
saccE = sacc[which(sacc$Type=="Exon"),]
saccJ = sacc[which(sacc$Type=="Jxn"),]
saccT = sacc[which(sacc$Type=="Tx"),]


pdf("sacc_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:12) {
	symi = saccG[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = saccG[i,"snps"]
	feati = saccG[i,"gene"]
	p_i = signif(saccG[i,"pvalue"],3)
	typei = saccG[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
for (i in 1:12) {
	symi = saccE[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = saccE[i,"snps"]
	feati = saccE[i,"gene"]
	p_i = signif(saccE[i,"pvalue"],3)
	typei = saccE[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
for (i in 1:12) {
	symi = saccJ[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = saccJ[i,"snps"]
	feati = saccJ[i,"gene"]
	p_i = signif(saccJ[i,"pvalue"],3)
	typei = saccJ[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
for (i in 1:12) {
	symi = saccT[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = saccT[i,"snps"]
	feati = saccT[i,"gene"]
	p_i = signif(saccT[i,"pvalue"],3)
	typei = saccT[i,"Type"]

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






## make CSV of top 1000 of each
sacc_merged = rbind(saccG[1:1000,],saccE[1:1000,],saccJ[1:1000,],saccT[1:1000,])
sacc_merged = sacc_merged[,-which(names(sacc_merged)=="gencodeTx")]

sacc = sacc_merged
sacc$EnsemblGeneID = ss(sacc$EnsemblGeneID, "\\.")

## snpMap
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
# snpMap = snpMap[which(rownames(snpMap) %in% c(sacc$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load("../../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../../data/zandiHypde_bipolar_rseGene_n511.rda")
gMap = as.data.frame(rowRanges(rse_gene))[,c("seqnames","start","end","strand","Class")]
eMap = as.data.frame(rowRanges(rse_exon))[,c("seqnames","start","end","strand","Class")]
jMap = as.data.frame(rowRanges(rse_jxn))[,c("seqnames","start","end","strand","Class")]
txMap = as.data.frame(rowRanges(rse_tx))[,c("seqnames","start","end","strand","source")]
txMap$source = "InGen"
# rm(rse_gene, rse_exon, rse_jxn, rse_tx)
colnames(gMap) = colnames(eMap) = colnames(jMap) = colnames(txMap) = 
	c("feat_chr","feat_start","feat_end","strand","Class")
featMap = rbind(rbind(rbind(gMap, eMap),jMap),txMap)
featMap$Type = c(rep("Gene",nrow(gMap)),rep("Exon",nrow(eMap)),rep("Jxn",nrow(jMap)),rep("Tx",nrow(txMap)))

geneMap = as.data.frame(rowRanges(rse_gene))[,c("gencodeID","Symbol","ensemblID","gene_type")]

## put together
snpMap_temp = snpMap[sacc$snps,]
featMap_temp = featMap[sacc$gene,]
geneMap_temp = geneMap[match(sacc$EnsemblGeneID, geneMap$ensemblID),]
sacc2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),sacc)
# sacc2 = sacc2[,-which(colnames(sacc2)=="gencodeTx")]

sacc3 = sacc2[,c(2,12:14,26,20,15:19,22:24,27:30)]

write.csv(sacc3, file="genomewide_snps_sacc_eqtls_top1000.csv")



