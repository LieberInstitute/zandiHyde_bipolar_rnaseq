######## 
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)

## load data
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
geneExprs = recount::getRPKM(rse_gene, "Length")

## deconvolution
load("/users/ajaffe/Lieber/Projects/RNAseq/MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")

yExprs_Scaled = scale(geneExprs[rownames(coefEsts),])

propEsts= minfi:::projectCellType(yExprs_Scaled, coefEsts)
propEstsScaled = prop.table(propEsts,1)
propEsts= as.data.frame(propEsts)

rse_gene$PrimaryDx[rse_gene$PrimaryDx == "Other"] = "Bipolar"
rse_gene$PrimaryDx = droplevels(rse_gene$PrimaryDx)

boxplot(propEsts$Neurons ~ rse_gene$PrimaryDx*rse_gene$BrainRegion,
	ylab = "Adult Neuron RNA Fraction")

boxplot(propEsts$Oligo ~ rse_gene$PrimaryDx*rse_gene$BrainRegion,
	ylab = "Adult Oligo RNA Fraction")
boxplot(propEsts$OPC ~ rse_gene$PrimaryDx*rse_gene$BrainRegion,
	ylab = "Adult OPC RNA Fraction")
boxplot(propEsts$OPC ~ rse_gene$PrimaryDx*rse_gene$BrainRegion,
	ylab = "Adult OPC RNA Fraction")

