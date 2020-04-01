##

library(MatrixEQTL)
library(GenomicRanges)
library(jaffelab)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")

# filter for age
aIndex= which(pd$Age > 13 & pd$Dx == "Control")
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex] )
jRpkm2 = as.matrix(log2(jRpkm[,aIndex]+1))
exonRpkm2= as.matrix(log2(exonRpkm[,aIndex]+1))
geneRpkm2= as.matrix(log2(geneRpkm[,aIndex]+1))
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

### filter rp80m
expIndex=which(rowMeans(jRpkm2) > 0.2 & jMap$code != "Novel")
jRpkm2 = jRpkm2[expIndex,]
jMap = jMap[expIndex]
rownames(snpMap) = rownames(snp2)
### load eQTLs
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/control/rdas/junction_eqtl_control_13plus_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/control/rdas/gene_eqtl_control_13plus_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/control/rdas/exon_eqtl_control_13plus_cisOnly.rda")

## which genes and snps

s = c("rs2075677", "rs4904738")
snpSub = snp2[grep(paste(s, collapse="|"),rownames(snp2)),]
snpMapSub = snpMap[grep(paste(s, collapse="|"),rownames(snpMap)),]

g = c("STAU1","LRFN5")
jIndex = which(jMap$newGeneSymbol %in% g)
jRpkmSub = jRpkm2[jIndex,]
jMapSub = jMap[jIndex,]
gIndex = which(geneMap$Symbol %in% g)
geneRpkmSub = geneRpkm2[gIndex,]
geneMapSub = geneMap[gIndex,]
eIndex = which(exonMap$Symbol %in% g)
exonRpkmSub = exonRpkm2[eIndex,]
exonMapSub = exonMap[eIndex,]

## clean
yExprsClean = rbind(cleaningY(geneRpkmSub,cbind(mod, pcsGene), P=1),
	cleaningY(exonRpkmSub, cbind(mod, pcsExon), P=1 ),
	cleaningY(jRpkmSub, cbind(mod, pcsJxn), P=1))

## rs4904738
rs4904738 = snpSub[1,]
rs4904738 = as.numeric(snpSub[1,])
rs4904738 = gsub(0, paste0(snpMapSub$ALT[1],
	snpMapSub$ALT[1]), rs4904738)
rs4904738 = gsub(1, paste0(snpMapSub$ALT[1], 
	snpMapSub$COUNTED[1]), rs4904738)
rs4904738 = gsub(2, paste0(snpMapSub$COUNTED[1], 
	snpMapSub$COUNTED[1]), rs4904738)
rs4904738 = factor(rs4904738, levels = c("TT", "CT", "CC")) # risk

## rs2075677
rs2075677 = snpSub[2,]
rs2075677 = as.numeric(snpSub[2,])
rs2075677 = gsub(0, paste0(snpMapSub$ALT[2],
	snpMapSub$ALT[2]), rs2075677)
rs2075677 = gsub(1, paste0(snpMapSub$ALT[2], 
	snpMapSub$COUNTED[2]), rs2075677)
rs2075677 = gsub(2, paste0(snpMapSub$COUNTED[2], 
	snpMapSub$COUNTED[2]), rs2075677)
rs2075677 = factor(rs2075677, levels = c("GG", "AG", "AA")) # risk

pdf("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/maynard/example_boxplots_rs4904738.pdf",useDingbats=FALSE)
par(mar=c(3,6,6,2), cex.axis=2, cex.lab=2, cex.main=1.7)
for(i in 1:nrow(yExprsClean)) {
	boxplot(yExprsClean[i,] ~ rs4904738, outline=FALSE,
		ylim= quantile(yExprsClean[i,],c(0.005, 0.995)), 
		ylab= "Adjusted log2 Expression",
		xlab="", main = rownames(yExprsClean)[i])
	points(yExprsClean[i,] ~ jitter(as.numeric(rs4904738),amount=0.1),
		pch = 21, bg="grey",cex=1.4)
}
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/maynard/example_boxplots_rs2075677.pdf",useDingbats=FALSE)
par(mar=c(3,6,6,2), cex.axis=2, cex.lab=2, cex.main=1.7)
for(i in c(37,80)) {
	boxplot(yExprsClean[i,] ~ rs2075677, outline=FALSE,
		ylim= quantile(yExprsClean[i,],c(0.005, 0.995)), 
		ylab= "Adjusted log2 Expression",
		xlab="", main = rownames(yExprsClean)[i])
	points(yExprsClean[i,] ~ jitter(as.numeric(rs2075677),amount=0.1),
		pch = 21, bg="grey",cex=1.4)
}
dev.off()
