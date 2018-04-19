##

library(jaffelab)
library(SummarizedExperiment)
library(readxl)
library(RColorBrewer)
# library(LIBDpheno)

## load phenotype data
pd = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/read_and_alignment_metrics_zandiHyde_Bipolar_LIBD.csv", stringsAsFactors=FALSE)
pd$FlowCell= ss(pd$SAMPLE_ID, "_",2)
colnames(pd)[55]= "NumTotalMapped"

aL=min(pd$totalAssignedGene)
aH=max(pd$totalAssignedGene)
mL=min(pd$overallMapRate)
mH=max(pd$overallMapRate)
mitoL=min(pd$mitoRate)
mitoH=max(pd$mitoRate)

# drop samples?
pdf("RIN_check_predrop.pdf", h=10,w=10)
par(mfcol=c(2,2),mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
plot(pd$RIN, pd$totalAssignedGene, pch =21, bg="grey")
abline(h=0.3, lty=2)
plot(pd$RIN, pd$mitoRate, pch =21, bg="grey")
plot(pd$RIN, pd$overallMapRate, pch =21, bg="grey")
abline(h=0.5, lty=2)
plot(pd$mitoRate, pd$overallMapRate, pch =21, bg="grey")
abline(h=0.5, lty=2)
abline(v=0.1, lty=2)
dev.off()

## drop 5 samples
pd = pd[pd$overallMapRate > 0.5 & pd$mitoRate < .1 & pd$totalAssignedGene > .3,]

pdf("RIN_check_postdrop.pdf", h=10,w=10)
par(mfcol=c(2,2),mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
plot(pd$RIN, pd$totalAssignedGene, pch =21, bg="grey", ylim=c(aL,aH))
abline(h=0.3, lty=2)
plot(pd$RIN, pd$mitoRate, pch =21, bg="grey", ylim=c(mitoL,mitoH))
plot(pd$RIN, pd$overallMapRate, pch =21, bg="grey", ylim=c(mL,mH) )
abline(h=0.5, lty=2)
plot(pd$mitoRate, pd$overallMapRate, pch =21, bg="grey", xlim=c(mitoL,mitoH), ylim=c(mL,mH))
abline(h=0.5, lty=2)
abline(v=0.1, lty=2)
dev.off()


##########
## fix some IDs based on genotyping

## drop 4 based on previous genotypes
dropIndex = which(pd$BrNum %in% c("Br1142","Br1143", 
	"Br1178","Br1179", "Br2407","Br2260", "Br2385","Br2538"))
pd = pd[-dropIndex,]



## load
load("/dcl01/lieber/ajaffe/lab/brain_swap/genotype_match_matrix_allBrains.rda")
snpCorMat = snpCor2[,match(pd$RNum, ss(colnames(snpCor2),"_"))]
colnames(snpCorMat) = pd$RNum

matchInd = as.data.frame(which(snpCorMat > 0.6, arr.ind=TRUE, useNames=TRUE))
matchInd$GenoInfo = rownames(snpCorMat)[matchInd$row]
matchInd$GenoBrNum = ss(matchInd$GenoInfo, "_")
matchInd$RnaInfo = colnames(snpCorMat)[matchInd$col]
matchInd$RnaBrNum = pd$BrNum[match(matchInd$RnaInfo, pd$RNum)]

## add correlations
for(i in 1:nrow(matchInd)) matchInd$genoCor[i] = snpCorMat[matchInd$row[i], matchInd$col[i]]
matchInd = matchInd[,c(7, 3:6)]

# check duplicates
matchInd[matchInd$RnaInfo %in% matchInd$RnaInfo[
	matchInd$RnaBrNum != matchInd$GenoBrNum],]
                    # genoCor        GenoInfo GenoBrNum RnaInfo RnaBrNum
# Br2473_Omni5M     0.8668090   Br2473_Omni5M    Br2473  R14071   Br2260
# Br2473_Omni5M.1   0.9033342   Br2473_Omni5M    Br2473  R14290   Br2260
# Br2260_Omni5M     0.8735736   Br2260_Omni5M    Br2260  R14077   Br2473
# Br2260_Omni5M.1   0.9029512   Br2260_Omni5M    Br2260  R14296   Br2473

# Br2301_Omni5M     0.9121581   Br2301_Omni5M    Br2301  R14996   Br2538
# Br2301_Omni5M.1   0.9127415   Br2301_Omni5M    Br2301  R15029   Br2538
# Br2538_Omni5M     0.9223550   Br2538_Omni5M    Br2538  R14103   Br2301

# Br2385_Omni5M     0.9237195   Br2385_Omni5M    Br2385  R14056   Br2533
# Br2385_Omni5M.1   0.9321742   Br2385_Omni5M    Br2385  R14992   Br2533

# Br5435_Macrogen   0.9637900 Br5435_Macrogen    Br5435  R14080   Br5434
# Br5434_Macrogen   0.9773414 Br5434_Macrogen    Br5434  R14087   Br5435
# Br5434_Macrogen.1 0.9215299 Br5434_Macrogen    Br5434  R14306   Br5435

### mismatch because of RNA mislabeling
# Br1148_h650       0.9766927     Br1148_h650    Br1148  R14210   Br1936
# Br5955_Macrogen   0.9385395 Br5955_Macrogen    Br5955  R15086   Br5944
# Br5962_Macrogen   0.8454689 Br5962_Macrogen    Br5962  R15087   Br5955
# Br5974_Macrogen   0.9904247 Br5974_Macrogen    Br5974  R15088   Br5962
# Br5979_Macrogen   0.9080365 Br5979_Macrogen    Br5979  R15089   Br5974









#### plots for Zandi 6/30/17

library(jaffelab)
library(SummarizedExperiment)
library(readxl)
library(RColorBrewer)
# library(LIBDpheno)

## ERCC data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rawCounts_zandiHyde_Bipolar_LIBD_n540.rda", verbose=TRUE)
rm(list=ls()[! ls() %in% c("erccTPM","metrics")])
pd=metrics
pd[,1:3] = lapply(pd[,1:3], function(x) as.character(x))
names(pd)[5] = "BrainRegion"
pd$SampleID = paste0(pd$BrNum, "_", pd$BrainRegion)
pd = pd[,c(1,70,2:69)]
pd$FlowCell= ss(pd$SAMPLE_ID, "_",2)
colnames(pd)[55]= "NumTotalMapped"



### plot
pdf("metrics_by_flowcell.pdf", h=5,w=5)
par(mar=c(7,6,2,2),cex.axis=1,cex.lab=1.5,cex.main=2)
palette(brewer.pal(8,"Dark2"))
boxplot(pd$overallMapRate ~ pd$FlowCell, las= 3,
          ylim = c(0,1),
          outline=FALSE,ylab="Alignment Rate")
points(pd$overallMapRate ~ jitter(as.numeric(
    factor(pd$FlowCell)),
	amount=0.15), pch=21,
    bg = as.numeric(factor(pd$BrainRegion)))
legend("bottomleft", levels(factor(pd$BrainRegion)),
       pch = 15, col = 1:2,cex=.8)
abline(h=.5, lty=2, col="grey")

boxplot(pd$totalAssignedGene ~ pd$FlowCell, las= 3,
          ylim = c(0,.7),
          outline=FALSE,ylab="Gene Assignement Rate")
points(pd$totalAssignedGene ~ jitter(as.numeric(
    factor(pd$FlowCell)),
	amount=0.15), pch=21,
    bg = as.numeric(factor(pd$BrainRegion)))
abline(h=.3, lty=2, col="grey")


dev.off()



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


#### pca 

## ercc pca
pcaERCC = prcomp(t(log2(erccTPM+1)))
pcaVarsERCC = getPcaVars(pcaERCC)


geneRpkm = getRPKM(rse_gene)
gMetrics = colData(rse_gene)
gMetrics$flowcell = ss(as.character(pd$SAMPLE_ID), "_",2)

geneRpkm = geneRpkm[which(rowMeans(geneRpkm) > .1),]
pca1 = prcomp(t(log2(geneRpkm+1)))
pcaVars1 = getPcaVars(pca1)

pdf("gene_pca.pdf", h=5,w=5)
par(mar=c(7,6,3,2),cex.axis=1,cex.lab=1.5,cex.main=1.5)
palette(brewer.pal(8,"Dark2"))
boxplot(pca1$x[,3] ~ gMetrics$flowcell, las= 3, main="Gene PCs",
          outline=FALSE,
		  ylab=paste0("PC3: ", pcaVars1[3], "% Var Expl"),
		  ylim=c(min(pca1$x[,3]),max(pca1$x[,3])))
points(pca1$x[,3] ~ jitter(as.numeric(
    factor(gMetrics$flowcell)),
	amount=0.15), pch=21,
    bg = as.numeric(factor(gMetrics$Brain.Region)))

boxplot(pcaERCC$x[,1] ~ pd$FlowCell, las= 3, main="ERCC PCs",
          outline=FALSE,
		  ylab=paste0("PC1: ", pcaVars1[1], "% Var Expl"))
points(pcaERCC$x[,1] ~ jitter(as.numeric(
    factor(pd$FlowCell)),
	amount=0.15), pch=21,
    bg = as.numeric(factor(pd$BrainRegion)))
dev.off()



### Regional differences ( in controls )
library(affy)

y = data.frame(Amyg=rowMeans(geneRpkm[,which(pd$BrainRegion=="Amygdala" & pd$PrimaryDx=="Control")]), 
				sacc=rowMeans(geneRpkm[,which(pd$BrainRegion=="sACC" & pd$PrimaryDx=="Control")]))
pdf("maplot.pdf")
ma.plot( rowMeans(log2(y)), log2(y[, 1])-log2(y[, 2]), cex=1 )
dev.off()




### we will switch some BrNums in the genotype data
##  swap Br2473 Br2260
##  swap Br2301 Br2538
##  
##  Br2385 --> Br2533
##  Br2538 --> Br2301
##  Br2407 --> Br1563

## switch RNA pheno info:
toxicant[[1]]$primarydx = as.character(toxicant[[1]]$primarydx)
toxicant[[1]]$race = as.character(toxicant[[1]]$race)
toxicant[[1]]$sex = as.character(toxicant[[1]]$sex)

swap = toxicant[[1]][toxicant[[1]]$brnumerical %in% c(1271,1281),
	c("brnum", "agedeath", "sex", "race", "primarydx", "ph","pmi")]

pd[pd$BrNum == "Br1721",c("BrNum", "age","Sex","Race","Dx","pH","PMI")] = rbind(swap[1,], swap[1,])
pd[pd$BrNum == "Br1289",c("BrNum", "age","Sex","Race","Dx","pH","PMI")] = rbind(swap[2,], swap[2,])

### check for whether matching
pdDrop = pd[!pd$RNum %in% matchInd$RnaInfo,] # one sample
# > unique(pdDrop$BrNum)
# [1] "Br2492" "Br5485" "Br5209" "Br1623" "Br1060" "Br1591"

pd = pd[pd$RNum %in% matchInd$RnaInfo,] # six brains

write.table(matchInd, file="qcChecks/genotype_checking_table.txt",
	sep="\t", row.names=FALSE, col.names=FALSE)

########################
### load counts ########

load("/dcl01/lieber/ajaffe/lab/brainseq_phase3/processed_data/rawCounts_Caudate_BrainSeq_Phase3_LIBD_n969.rda")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase3/processed_data/rpkmCounts_Caudate_BrainSeq_Phase3_LIBD_n969.rda")
rm(jRpkm)

## filter to same samples
geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]
# filter

jCounts = jCounts[,pd$SAMPLE_ID]

# 1 or more reads in 10 or more samples
jCountsLogical = DataFrame(sapply(jCounts, function(x) x > 0))
jIndex = which(rowSums(as.data.frame(jCountsLogical)) >= 10) 

jCounts =  as.matrix(as.data.frame(jCounts[jIndex,]))
jMap = jMap[jIndex,]

## update junction annotation class
tt = rowSums(!is.na(as.data.frame(mcols(jMap)[,
	c("startExon", "endExon")])))
jMap$Class[tt == 1] = "AltStartEnd"
jMap$Class[tt == 2 & !jMap$inGencode] = "ExonSkip"

txTpm = txTpm[,pd$SAMPLE_ID]

##################
## compare full counts
pca = prcomp(t(log2(geneCounts+1)))
plot(pca$x) # looks fine

###############
### merge
sIndexes = splitit(pd$RNum)
geneCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(geneCounts[,ii])))
})

exonCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(exonCounts[,ii])))
})

jCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(jCounts[,ii])))
})

## transcripts
combineTxQuant = function (sampID1, sampID2, salmonDir)
{    
    ##get names of transcripts
    txNames = read.table(file.path(salmonDir, sampID1, "quant.sf"),header = TRUE)$Name
    txNames = as.character(txNames)
    gencodeTx = ss(txNames, "\\|",1)
    
    txQuant1 = read.table(file.path(salmonDir, sampID1, "quant.sf"),header = TRUE)
    txQuant2 = read.table(file.path(salmonDir, sampID2, "quant.sf"),header = TRUE)
    
    txQuant = data.frame(EffectiveLength=rep(NA,length(gencodeTx)), NumReads=NA, RPK=NA, TPM=NA, row.names=gencodeTx)
    
    txQuant$EffectiveLength = txQuant1$EffectiveLength + txQuant2$EffectiveLength
    txQuant$NumReads = txQuant1$NumReads + txQuant2$NumReads
    txQuant$RPK = txQuant$NumReads/txQuant$EffectiveLength
    txQuant$TPM = 1e6 * (txQuant$NumReads/txQuant$EffectiveLength) / sum(txQuant$RPK)

    return(txQuant)
}

doubleSamples = lapply(sIndexes[lengths(sIndexes)>1], function(ii) pd$SAMPLE_ID[ii])
doubleSamples = do.call("rbind",doubleSamples)

## add transcripts
txPath = "/dcl01/lieber/ajaffe/lab/brainseq_phase3/processed_data/Salmon_tx/"
doubleTxCombine = apply(doubleSamples, 1, function(x) {
	cat(".")
	combineTxQuant(x[1], x[2], txPath)
})
doubleTx = sapply(doubleTxCombine, "[[", "TPM")

singleTx = txTpm[,unlist(sIndexes[lengths(sIndexes) == 1])]
tpmMatMerge = cbind(singleTx, doubleTx)
colnames(tpmMatMerge) = ss(colnames(tpmMatMerge), "_")
tpmMatMerge = tpmMatMerge[,names(sIndexes)]

############
## make expression sets
library(SummarizedExperiment)
pdDF = DataFrame(pd)
classIndex = splitit(sapply(pdDF, class))

pdMergeList = sapply(sIndexes, function(ii) {
	x = pdDF[ii,]
	o = pdDF[ii[1],]
	for(j in classIndex$character) o[j] = CharacterList(unique(x[,j]))
	for(j in classIndex$integer) o[j] = IntegerList(unique(x[,j]))
	for(j in classIndex$logical) o[j] = LogicalList(unique(x[,j]))
	for(j in classIndex$numeric) o[j] = NumericList(unique(x[,j]))
	return(o)
})
pdMerge = do.call("rbind", pdMergeList)		

## fix columns
colnames(pdMerge)[5] = "Age"
pdMerge$RNum = sapply(pdMerge$RNum, "[[", 1)
pdMerge$BrNum = sapply(pdMerge$BrNum, "[[", 1)
pdMerge$Region = sapply(pdMerge$Region, "[[", 1)
pdMerge$RIN = sapply(pdMerge$RIN, "[[", 1)
pdMerge$Age = sapply(pdMerge$Age, "[[", 1)
pdMerge$Sex = sapply(pdMerge$Sex, "[[", 1)
pdMerge$Race = sapply(pdMerge$Race, "[[", 1)
pdMerge$Dx = sapply(pdMerge$Dx, "[[", 1)
pdMerge$pH = sapply(pdMerge$pH, "[[", 1)
pdMerge$PMI = sapply(pdMerge$PMI, "[[", 1)

rownames(pdMerge) = pdMerge$RNum
identical(rownames(pdMerge), colnames(geneCountMerge))

##########
## gene
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$gencodeTx = CharacterList(strsplit(geneMapGR$gencodeTx, ";"))
rse_gene = SummarizedExperiment(
	assays = list('counts' = geneCountMerge),
    colData = pdMerge, rowRanges = geneMapGR)

## gene
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
exonMapGR$gencodeTx = CharacterList(strsplit(exonMapGR$gencodeTx, ";"))
rse_exon = SummarizedExperiment(
	assays = list('counts' = exonCountMerge),
    colData = pdMerge, rowRanges = exonMapGR)
	
## junction
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCountMerge),
    colData = pdMerge, rowRanges = jMap)
	
## transcript
gtf = rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
tx = gtf[which(gtf$type == "transcript")]
names(tx) = tx$transcript_id
txMap = tx[rownames(tpmMatMerge)]
rse_tx = SummarizedExperiment(
	assays = list('tpm' = tpmMatMerge),
    colData = pdMerge, rowRanges = txMap)
	
#### function for RPKM
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}

getRPM = function(rse, target = 80e6) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped/target), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

################
## save ########
dir.create("count_data")
save(rse_gene, getRPKM, compress=TRUE,
	file = "count_data/caudate_brainseq_phase3_hg38_rseGene_merged_n466.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "count_data/caudate_brainseq_phase3_hg38_rseExon_merged_n466.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "count_data/caudate_brainseq_phase3_hg38_rseJxn_merged_n466.rda")
save(rse_tx, compress=TRUE,
	file = "count_data/caudate_brainseq_phase3_hg38_rseTx_merged_n466.rda")

