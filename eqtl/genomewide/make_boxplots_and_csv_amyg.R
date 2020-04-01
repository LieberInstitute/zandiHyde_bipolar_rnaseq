####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")

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
load("mergedEqtl_output_amyg_genomewide_4features_FDR01.rda", verbose=TRUE)
amyg = allEqtlFDR01


##
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))
mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)


################
## load expression
## amygdala
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
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
amyg$Symbol = as.character(amyg$Symbol)

amygG = amyg[which(amyg$Type=="Gene"),]
amygE = amyg[which(amyg$Type=="Exon"),]
amygJ = amyg[which(amyg$Type=="Jxn"),]
amygT = amyg[which(amyg$Type=="Tx"),]


pdf("amyg_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:12) {
	symi = amygG[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = amygG[i,"snps"]
	feati = amygG[i,"gene"]
	p_i = signif(amygG[i,"pvalue"],3)
	typei = amygG[i,"Type"]

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
	symi = amygE[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = amygE[i,"snps"]
	feati = amygE[i,"gene"]
	p_i = signif(amygE[i,"pvalue"],3)
	typei = amygE[i,"Type"]

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
	symi = amygJ[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = amygJ[i,"snps"]
	feati = amygJ[i,"gene"]
	p_i = signif(amygJ[i,"pvalue"],3)
	typei = amygJ[i,"Type"]

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
	symi = amygT[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = amygT[i,"snps"]
	feati = amygT[i,"gene"]
	p_i = signif(amygT[i,"pvalue"],3)
	typei = amygT[i,"Type"]

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
amyg_merged = rbind(amygG[1:1000,],amygE[1:1000,],amygJ[1:1000,],amygT[1:1000,])
amyg_merged = amyg_merged[,-which(names(amyg_merged)=="gencodeTx")]

amyg = amyg_merged
amyg$EnsemblGeneID = ss(amyg$EnsemblGeneID, "\\.")

## snpMap
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
# snpMap = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

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
snpMap_temp = snpMap[amyg$snps,]
featMap_temp = featMap[amyg$gene,]
geneMap_temp = geneMap[match(amyg$EnsemblGeneID, geneMap$ensemblID),]
amyg2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),amyg)
# amyg2 = amyg2[,-which(colnames(amyg2)=="gencodeTx")]

amyg3 = amyg2[,c(2,12:14,26,20,15:19,22:24,27:30)]

write.csv(amyg3, file="genomewide_snps_amyg_eqtls_top1000.csv")



