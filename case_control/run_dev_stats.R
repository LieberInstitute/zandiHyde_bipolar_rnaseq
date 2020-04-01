####
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(recount)
library(jaffelab)
library(IRanges)

## feature maps
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseTx_n511.rda")


## common expression map to annotate
rowRanges(rse_gene)$gencodeID = rownames(rse_gene)
rowRanges(rse_gene)$Class = "InGen"
rowRanges(rse_gene)$Type = "Gene"

rowRanges(rse_exon)$Class = "InGen"
rowRanges(rse_exon)$Type = "Exon"

colnames(mcols(rowRanges(rse_jxn)))[13:14] = c("gencodeID", "Symbol")
rowRanges(rse_jxn)$Type = "Junction"

colnames(mcols(rowRanges(rse_tx)))[c(5,8)] = c("gencodeID", "Symbol")
rowRanges(rse_tx)$Class = "InGen"
rowRanges(rse_tx)$Type = "Transcript"

name = c("gencodeID", "Symbol", "Type", "Class")
exprsMap = rbind(as.data.frame(rowRanges(rse_gene))[,name],
	as.data.frame(rowRanges(rse_exon))[,name],
	as.data.frame(rowRanges(rse_jxn))[,name],
	as.data.frame(rowRanges(rse_tx))[,name])
exprsMap = DataFrame(exprsMap)

## add entrez
exprsMap$EntrezID = rowData(rse_gene)$EntrezID[
	match(exprsMap$gencodeID, rownames(rse_gene))]

################
## load

load("bipolarControl_deStats_byRegion_qSVAjoint.rda", verbose=TRUE)
exprsMap = exprsMap[rownames(statOut),]
statOut = cbind(statOut, exprsMap)


###############################
## add dev stats from phase2 ##

devObjs = list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda",
	pattern = "_only", full = TRUE)
devObjs = devObjs[c(3,1,5,7,4,2,6,8)]
names(devObjs) = c("gene_DLPFC","exon_DLPFC","jxn_DLPFC","tx_DLPFC",
	"gene_HIPPO","exon_HIPPO","jxn_HIPPO","tx_HIPPO")
devStats = mclapply(devObjs, function(x) {
	cat(".")
	load(x)
	top$Feature = rownames(top)
	return(top)
},mc.cores = 4)

## add bonf
devStats = lapply(devStats, function(x) {
	x$p_bonf = p.adjust(x$P.Value, method = "bonf")
	return(x)
})
sapply(devStats, function(x) mean(x$p_bonf < 0.05))

## merge
devStatsDLPFC = do.call("rbind", devStats[grep("DLPFC", names(devStats))])
devStatsHIPPO = do.call("rbind", devStats[grep("HIPPO", names(devStats))])
rownames(devStatsDLPFC) = devStatsDLPFC$Feature
rownames(devStatsHIPPO) = devStatsHIPPO$Feature
identical(rownames(devStatsHIPPO), rownames(devStatsDLPFC))

## combine with DE
devStatsDLPFC = devStatsDLPFC[match(rownames(statOut), rownames(devStatsDLPFC)),]
devStatsHIPPO = devStatsHIPPO[match(rownames(statOut), rownames(devStatsHIPPO)),]

statOut$inDevData = !is.na(devStatsDLPFC$F)
statOut$devReg_DLPFC = devStatsDLPFC$adj.P.Val < 0.05
statOut$devReg_HIPPO = devStatsHIPPO$adj.P.Val < 0.05
statOut$negCorr_DLPFC = devStatsDLPFC$adj.P.Val < 0.05 & devStatsDLPFC$ageCorr < 0
statOut$negCorr_HIPPO = devStatsHIPPO$adj.P.Val < 0.05 & devStatsHIPPO$ageCorr < 0
statOut$posCorr_DLPFC = devStatsDLPFC$adj.P.Val < 0.05 & devStatsDLPFC$ageCorr > 0
statOut$posCorr_HIPPO = devStatsHIPPO$adj.P.Val < 0.05 & devStatsHIPPO$ageCorr > 0

## splitup
statList = split(statOut, statOut$Type)

##################################
### find overlaps by features ####
##################################

## amygdala
oStats_Amyg = lapply(statList, function(x) {
	o = rbind(summary(lm(x$t_Amyg ~ x$negCorr_DLPFC))$coef[2,c(1,4)],
		summary(lm(x$t_Amyg ~ x$posCorr_DLPFC))$coef[2,c(1,4)],
		summary(lm(x$t_Amyg ~ x$negCorr_HIPPO))$coef[2,c(1,4)],
		summary(lm(x$t_Amyg ~ x$posCorr_HIPPO))$coef[2,c(1,4)])
	rownames(o) = c("fetal_DLPFC","postnatal_DLPFC",
		"fetal_HIPPO","postnatal_HIPPO")
	return(o)
})

## sACC
oStats_sACC = lapply(statList, function(x) {
	o = rbind(summary(lm(x$t_sACC ~ x$negCorr_DLPFC))$coef[2,c(1,4)],
		summary(lm(x$t_sACC ~ x$posCorr_DLPFC))$coef[2,c(1,4)],
		summary(lm(x$t_sACC ~ x$negCorr_HIPPO))$coef[2,c(1,4)],
		summary(lm(x$t_sACC ~ x$posCorr_HIPPO))$coef[2,c(1,4)])
	rownames(o) = c("fetal_DLPFC","postnatal_DLPFC",
		"fetal_HIPPO","postnatal_HIPPO")
	return(o)
})

oStats_Amyg[c("Gene", "Exon", "Transcript", "Junction")]
oStats_sACC[c("Gene", "Exon", "Transcript", "Junction")]
