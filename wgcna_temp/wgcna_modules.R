##
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)

## load data
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/degradation_rse_BipSeq_BothRegions.rda")

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
			
## add ancestry 
load("../genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

###########
# filter ##
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

geneMap = as.data.frame(rowRanges(rse_gene))


#########################
## load wgcna output ####
#########################

## load
load("rdas/constructed_network_signed_bicor.rda")

# get colors
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]


gList = split(rowData(rse_gene)$gencodeID, 
	factor(net$colorsLab,levels = colorDat$col))

	
red = data.frame(gencodeID=gList[["red"]])
pink = data.frame(gencodeID=gList[["pink"]])
magenta = data.frame(gencodeID=gList[["magenta"]])
royalblue = data.frame(gencodeID=gList[["royalblue"]])

red = geneMap[as.character(red$gencodeID),1:12]
pink = geneMap[as.character(pink$gencodeID),1:12]
magenta = geneMap[as.character(magenta$gencodeID),1:12]
royalblue = geneMap[as.character(royalblue$gencodeID),1:12]

write.csv(red, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna_temp/wgcna_red.csv")
write.csv(pink, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna_temp/wgcna_pink.csv")
write.csv(magenta, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna_temp/wgcna_magenta.csv")
write.csv(royalblue, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna_temp/wgcna_royalblue.csv")

	