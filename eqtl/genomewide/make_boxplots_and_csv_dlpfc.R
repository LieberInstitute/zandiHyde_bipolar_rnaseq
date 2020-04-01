####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pd = colData(rse_gene)

## load SNP data
load("rdas/overlappingSNPs.rda")  # snpMapKeep
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda")
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
load("mergedEqtl_output_dlpfc_genomewide_4features_FDR01.rda", verbose=TRUE)
dlp = allEqtlFDR01


##
pd$Dx = factor(pd$Dx,
	levels = c("Control", "Bipolar"))
mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)


################
## load expression
## dlpfc
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
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
dlp$Symbol = as.character(dlp$Symbol)

dlpG = dlp[which(dlp$Type=="Gene"),]
dlpE = dlp[which(dlp$Type=="Exon"),]
dlpJ = dlp[which(dlp$Type=="Jxn"),]
dlpT = dlp[which(dlp$Type=="Tx"),]


pdf("dlpfc_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:12) {
	symi = dlpG[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlpG[i,"snps"]
	feati = dlpG[i,"gene"]
	p_i = signif(dlpG[i,"pvalue"],3)
	typei = dlpG[i,"Type"]

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
	symi = dlpE[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlpE[i,"snps"]
	feati = dlpE[i,"gene"]
	p_i = signif(dlpE[i,"pvalue"],3)
	typei = dlpE[i,"Type"]

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
	symi = dlpJ[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlpJ[i,"snps"]
	feati = dlpJ[i,"gene"]
	p_i = signif(dlpJ[i,"pvalue"],3)
	typei = dlpJ[i,"Type"]

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
	symi = dlpT[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlpT[i,"snps"]
	feati = dlpT[i,"gene"]
	p_i = signif(dlpT[i,"pvalue"],3)
	typei = dlpT[i,"Type"]

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
dlp_merged = rbind(dlpG[1:1000,],dlpE[1:1000,],dlpJ[1:1000,],dlpT[1:1000,])
dlp_merged = dlp_merged[,-which(names(dlp_merged)=="gencodeTx")]

dlp = dlp_merged
dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")

## load SNP data
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda", verbose=TRUE)
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
# snpMap2 = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
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


snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
# dlp2 = dlp2[,-which(colnames(dlp2)=="gencodeTx")]

dlp3 = dlp2[,c(2,12,13,14,26,20,15:19,22:24,27:30)]
write.csv(dlp3, file="genomewide_snps_dlpfc_eqtls_top1000.csv")



