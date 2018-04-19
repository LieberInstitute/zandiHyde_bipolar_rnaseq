## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)


###############################################################
######################## load genomic data ####################

### transcript TPM
load("../preprocessed_data/rpkmCounts_zandiHyde_Bipolar_LIBD_n540.rda", verbose=TRUE)
rm(list=ls()[! ls() %in% c("txTpm","txMap")])

### counts
load("../preprocessed_data/rawCounts_zandiHyde_Bipolar_LIBD_n540.rda", verbose=TRUE)

## clean pd
pd = metrics
pd[,1:3] = lapply(pd[,1:3], function(x) as.character(x))
names(pd)[5] = "BrainRegion"
pd$SampleID = paste0(pd$BrNum, "_", pd$BrainRegion)
pd = pd[,c(1,70,2:69)]
##

jCounts = as.matrix(as.data.frame(jCounts))



###############################################################
######################## drop samples #########################

#####################
## drop based on alignment metrics  - drops 5
keepIndex = which(pd$overallMapRate>0.5 & pd$mitoRate<.1 & pd$totalAssignedGene>.3)
pd = pd[keepIndex,]

#####################
## drop based on inconsistent genotype data  - drops 8
## note Br1697_Amygdala already dropped in previous step (oMapRate=0.32)
droplist1 = c("Br1936_sACC","Br5974_Amygdala","Br1443_sACC","Br5901_Amygdala","Br1697_Amygdala",
			"Br5168_sACC","Br5205_Amygdala","Br5266_sACC","Br5434_sACC")
dropIndex1 = which(pd$SampleID %in% droplist1)
pd = pd[-dropIndex1,]

#####################
### drop based on region identity check  - drops 14
### see /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/qc_checks/check_regions.R
droplist2 = c("R14081","R14016","R15039","R15115","R13925","R14179","R14295",
			"R14228","R14175","R14269","R14235","R14231","R14151","R14242")			
dropIndex2 = which(pd$RNum %in% droplist2)
pd = pd[-dropIndex2,]

#####################
### drop brain Br1836 - low quality array per Andrew
droplist3 = "Br1836"		
dropIndex3 = which(pd$BrNum==droplist3)
pd = pd[-dropIndex3,]


geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]
jCounts = jCounts[,pd$SAMPLE_ID]





###############################################################
######################## update labels ########################

#####################
### update BrNumbers based on genotype data

# Br5962_Amygdala should be Br5974_Amygdala
ind = which(pd$SampleID == "Br5962_Amygdala")
pd[ind,c(2,4)] = c("Br5974_Amygdala","Br5974")

# Br5955_Amygdala should be Br5962_Amygdala
ind = which(pd$SampleID == "Br5955_Amygdala")
pd[ind,c(2,4)] = c("Br5962_Amygdala","Br5962")

# Br5944_Amygdala should be Br5955_Amygdala
ind = which(pd$SampleID == "Br5944_Amygdala")
pd[ind,c(2,4)] = c("Br5955_Amygdala","Br5955")

# Rename Br0922 to Br922 to match genotype labels
ind = which(pd$SampleID=="Br0922_sACC")
pd[ind,c(2,4)] = c("Br922_sACC","Br922")


#####################
### update Regions based on region identity check
### see /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/qc_checks/check_regions.R

### swap Br1562 regions
ind = which(pd$BrNum == "Br1562")
pd[ind,2] = c("Br1562_sACC","Br1562_Amygdala")
pd[ind,6] = c("sACC","Amygdala")

### R15072 - Br5939_Amygdala should be Br5939_sACC 
ind = which(pd$RNum == "R15072")
pd[ind,2] = "Br5939_sACC"
pd[ind,6] = "sACC"




###############################################################
######################## save data ############################

save(pd, file="annotated_phenotype_data_zandi_hyde_bipolar_n511.rda")

############
## make expression sets
rownames(pd) = pd$SampleID
pdDF = DataFrame(pd)
identical(pdDF$SAMPLE_ID, colnames(geneCounts))

## gene
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$gencodeTx = CharacterList(strsplit(geneMapGR$gencodeTx, ";"))
colnames(geneCounts) = rownames(pdDF)
rse_gene = SummarizedExperiment(
	assays = list('counts' = geneCounts),
    colData = pdDF, rowRanges = geneMapGR)

## exon
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
exonMapGR$gencodeTx = CharacterList(strsplit(exonMapGR$gencodeTx, ";"))
colnames(exonCounts) = rownames(pdDF)
rse_exon = SummarizedExperiment(
	assays = list('counts' = exonCounts),
    colData = pdDF, rowRanges = exonMapGR)
	
## junction
jIndex = which(rowSums(jCounts > 0) > 5)  ## present in >1% of samples
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = rownames(pdDF)
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCounts),
    colData = pdDF, rowRanges = jMap)
	
## transcript
gtf = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
tx = gtf[which(gtf$type == "transcript")]
names(tx) = tx$transcript_id
identical(names(tx), rownames(txTpm))

txTpm = txTpm[,pd$SAMPLE_ID]    ## filter to n=511
colnames(txTpm) = rownames(pdDF)
rse_tx = SummarizedExperiment(
	assays = list('tpm' = txTpm),
    colData = pdDF, rowRanges = tx)
	
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

##### update mean expression
gRpkm = getRPKM(rse_gene)
geneMapGR$meanExprs = rowMeans(gRpkm)
rse_gene = SummarizedExperiment(
	assays = list('counts' = geneCounts),
    colData = pdDF, rowRanges = geneMapGR)

eRpkm = getRPKM(rse_exon)
exonMapGR$meanExprs = rowMeans(eRpkm)
rse_exon = SummarizedExperiment(
	assays = list('counts' = exonCounts),
    colData = pdDF, rowRanges = exonMapGR)
	
jRpkm = getRPM(rse_jxn)
jMap$meanExprs = rowMeans(jRpkm)
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCounts),
    colData = pdDF, rowRanges = jMap)

	
################
## save ########

save(rse_gene, getRPKM, compress=TRUE,
	file = "zandiHypde_bipolar_rseGene_n511.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "zandiHypde_bipolar_rseExon_n511.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "zandiHypde_bipolar_rseJxn_n511.rda")
save(rse_tx, compress=TRUE,
	file = "zandiHypde_bipolar_rseTx_n511.rda")

