################################## check brain regions
## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)
library(RColorBrewer)
#### function for RPKM
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata")
## clean pd
pd = colData(rse_gene)
pd[,1:3] = lapply(pd[,1:3], function(x) as.character(x))
names(pd)[5] = "BrainRegion"
pd$SampleID = paste0(pd$BrNum, "_", pd$BrainRegion)
pd = pd[,c(1,70,2:69)]

gRpkm = getRPKM(rse_gene)

#######
## drop 
keepIndex = which(pd$overallMapRate>0.5 & pd$mitoRate<.1 & pd$totalAssignedGene>.3)
pd = pd[keepIndex,]
gRpkm = gRpkm[,pd$SAMPLE_ID]

#######
## drop     # note Br1697_Amygdala already dropped in previous step (oMapRate=.32)
droplist = c("Br1936_sACC","Br5974_Amygdala","Br1443_sACC","Br5901_Amygdala","Br1697_Amygdala",
			"Br5168_sACC","Br5205_Amygdala","Br5266_sACC","Br5434_sACC")
dropIndex = which(pd$SampleID %in% droplist)
pd = pd[-dropIndex,]
gRpkm = gRpkm[,pd$SAMPLE_ID]


## filter low expressing genes
gRpkm = gRpkm[which(rowMeans(gRpkm)>0.1),]
yExprs = log2(gRpkm+1)

### top 1000 genes different between regions
library(genefilter)
p = rowttests(yExprs, pd$BrainRegion)
pOrd = p[order(p$p.value)[1:1000],]
ind1000 = which(rownames(p) %in% rownames(pOrd))
yExprs1000 = yExprs[ind1000,]

### estimate brain region
amyg = ifelse(pd$BrainRegion=="Amygdala",1,0)
sacc = ifelse(pd$BrainRegion=="sACC",1,0)
mod = data.frame(model.matrix(~amyg+sacc - 1))

library(limma)
fit = lmFit(yExprs1000, mod)
Xmat = fit$coef
Dmat = t(Xmat)%*%Xmat
guess = apply(yExprs1000, 2, function(x)  solve(Dmat, t(Xmat) %*% x))[2,]
identical(names(guess),pd$SAMPLE_ID)
pd$guess = guess

### plot
pdf("region_check_1000.pdf", h=6,w=5)
par(mar=c(8,6,4,2),cex.axis=1.5,cex.lab=1.5,cex.main=2)
palette(brewer.pal(8,"Dark2"))
boxplot(pd$guess ~ pd$BrainRegion, las= 3,
          ylim = range(pd$guess),
          outline=FALSE,ylab="sACC Identity")
points(pd$guess ~ jitter(as.numeric(
    factor(pd$BrainRegion)),
	amount=0.15), pch=21,
    bg = as.numeric(factor(pd$BrainRegion)))
segments(0,0.57,1.5,0.57, lty=2, col="grey")
segments(1.5,0.25,2.5,0.25, lty=2, col="grey")
dev.off()

pdf("region_check_3.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.4,cex.lab=1.4,cex.main=2)
palette(brewer.pal(8,"Dark2"))
plot(pd$totalAssignedGene, pd$guess, pch = 21, 
	 bg=as.numeric(factor(pd$BrainRegion)),
	 cex=1.5, main="",
     xlab="Gene Assignment Rate",
     ylab="")
legend("topleft", paste0(levels(factor(pd$BrainRegion))),
       pch = 15, col = 1:8,cex=.9)
abline(h=c(.3,.7), lty=2, col="grey")
dev.off()


pd1 = pd[which((pd$guess>0.57 & pd$BrainRegion=="Amygdala") | (pd$guess<0.25 & pd$BrainRegion=="sACC")),c(1:6,71)]
pd1 = as.data.frame(pd1[order(pd1$BrainRegion,pd1$guess),])
pd1
                        # SAMPLE_ID        SampleID   RNum  BrNum RIN BrainRegion        guess
# R13925_H7K5NBBXX R13925_H7K5NBBXX Br1458_Amygdala R13925 Br1458 7.0    Amygdala  0.610068006
# R15115_HFFGHBBXX R15115_HFFGHBBXX Br5799_Amygdala R15115 Br5799 8.6    Amygdala  0.650380929
# R15039_HFFGHBBXX R15039_HFFGHBBXX Br5826_Amygdala R15039 Br5826 6.5    Amygdala  0.743235822
# R13897_HCTYLBBXX R13897_HCTYLBBXX Br1562_Amygdala R13897 Br1562 8.0    Amygdala  0.790226619
# R14016_H7L3FBBXX R14016_H7L3FBBXX Br5165_Amygdala R14016 Br5165 6.7    Amygdala  0.872747866
# R14081_H7JM5BBXX R14081_H7JM5BBXX Br5611_Amygdala R14081 Br5611 7.0    Amygdala  0.941233185
# R15072_HFF3TBBXX R15072_HFF3TBBXX Br5939_Amygdala R15072 Br5939 7.5    Amygdala  1.052423270
# R14179_HCTYLBBXX R14179_HCTYLBBXX     Br1469_sACC R14179 Br1469 7.6        sACC -0.332851687
# R14295_H7JKMBBXX R14295_H7JKMBBXX     Br2071_sACC R14295 Br2071 6.5        sACC -0.329935229
# R14228_H7JKMBBXX R14228_H7JKMBBXX     Br2333_sACC R14228 Br2333 6.2        sACC -0.274313955
# R14175_H7JKMBBXX R14175_H7JKMBBXX     Br1752_sACC R14175 Br1752 7.5        sACC -0.190272580
# R14269_H7JKMBBXX R14269_H7JKMBBXX     Br5555_sACC R14269 Br5555 6.5        sACC -0.108019973
# R14235_HCTYLBBXX R14235_HCTYLBBXX     Br1661_sACC R14235 Br1661 6.2        sACC -0.027656924
# R14110_H7L3FBBXX R14110_H7L3FBBXX     Br1562_sACC R14110 Br1562 7.1        sACC -0.008095872
# R14231_H7JLCBBXX R14231_H7JLCBBXX     Br2454_sACC R14231 Br2454 6.7        sACC  0.093006560
# R14151_H7JHNBBXX R14151_H7JHNBBXX     Br2582_sACC R14151 Br2582 7.9        sACC  0.134227500
# R14242_HCTYLBBXX R14242_HCTYLBBXX     Br5190_sACC R14242 Br5190 5.9        sACC  0.198791829


### swap Br1562 regions
ind = which(pd$BrNum == "Br1562")
pd[ind,2] = c("Br1562_sACC","Br1562_Amygdala")
pd[ind,6] = c("sACC","Amygdala")

### R15072 - Br5939_Amygdala should be Br5939_sACC 
ind = which(pd$RNum == "R15072")
pd[ind,2] = "Br5939_sACC"
pd[ind,6] = "sACC"

### Drop other 14
ind1 = which(pd1$BrNum %in% c("Br5939","Br1562"))
pd2 = pd1[-ind1,]
dropInd = which(pd$RNum %in% pd2$RNum)
pd = pd[-dropInd,]







################################## check Br5939  -- PCA
## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)
library(RColorBrewer)
#### function for RPKM
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata")
## clean pd
pd = colData(rse_gene)
pd[,1:3] = lapply(pd[,1:3], function(x) as.character(x))
names(pd)[5] = "BrainRegion"
pd$SampleID = paste0(pd$BrNum, "_", pd$BrainRegion)
pd = pd[,c(1,70,2:69)]

gRpkm = getRPKM(rse_gene)

#######
## drop 
keepIndex = which(pd$overallMapRate>0.5 & pd$mitoRate<.1 & pd$totalAssignedGene>.3)
pd = pd[keepIndex,]
gRpkm = gRpkm[,pd$SAMPLE_ID]

#######
## drop 
droplist = c("Br1936_sACC","Br5974_Amygdala","Br1443_sACC","Br5901_Amygdala","Br1697_Amygdala",
			"Br5168_sACC","Br5205_Amygdala","Br5266_sACC","Br5434_sACC")
dropIndex = which(pd$SampleID %in% droplist)
pd = pd[-dropIndex,]
gRpkm = gRpkm[,pd$SAMPLE_ID]

## filter
gRpkm = gRpkm[which(rowMeans(gRpkm)>0.1),]

## get top 100 genes
library(genefilter)
p = rowttests(gRpkm, pd$BrainRegion)
pOrd = p[order(p$p.value)[1:100],]
ind100 = which(rownames(p) %in% rownames(pOrd))
gRpkm100 = gRpkm[ind100,]

## PCA
pca1 = prcomp(t(log2(gRpkm100+1)))
pcaVars1 = getPcaVars(pca1)

## PC 1 vs 2: Different brain regions
pdf("pca_log2Rpkm_PC1_2_100genes.pdf", h=6,w=12)
par(mfcol=c(1,2),mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"Dark2"))
plot(pca1$x[,1], pca1$x[,2], pch = as.numeric(pd$BrainRegion)+20, 
	 bg=ifelse(pd$BrNum=="Br5939","yellow",as.numeric(factor(pd$BrainRegion))),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("bottomleft", paste0(levels(factor(pd$BrainRegion))),
       pch = c(16,15), col = 1:8,cex=.9)
legend("bottomright", "Br5939", pch=22, col="gold", pt.bg="yellow", cex=.9)	  
plot(pca1$x[,1], pca1$x[,2], pch = as.numeric(pd$BrainRegion)+20, 
	 bg=ifelse(pd$BrNum=="Br5611","yellow",as.numeric(factor(pd$BrainRegion))),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("bottomleft", paste0(levels(factor(pd$BrainRegion))),
       pch = c(16,15), col = 1:8,cex=.9)	  
legend("bottomright", "Br5611", pch=22, col="gold",pt.bg="yellow", cex=.9)	     
dev.off()

pd1 = pd
pd1$PC1 = pca1$x[,1]
pd1[which(pd1$BrainRegion=="Amygdala" & pd1$PC1>1),c(1:8,71)]


