##### 
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)

## load data
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/degradation_rse_BipSeq_BothRegions.rda")
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

## make factor
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
				
## add ancestry 
load("../genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## gene filter
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]


#########################
## split back by region #
#########################

# model w/o interaction to subset by region
modSep = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

#### both ####
sACC_Index = which(colData(rse_gene)$BrainRegion == "sACC")
mod_sACC = modSep[sACC_Index,]
Amyg_Index = which(colData(rse_gene)$BrainRegion == "Amygdala")
mod_Amyg = modSep[Amyg_Index,]

#################
### Gene ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_gene[,sACC_Index])$counts, 
	genes = rowData(rse_gene))
dge_sACC = calcNormFactors(dge_sACC)
vGene_sACC = voom(dge_sACC,mod_sACC, plot=TRUE)

fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outGene_sACC = topTable(eBGene_sACC,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene_sACC = outGene_sACC[rownames(rse_gene),]
sum(outGene_sACC$adj.P.Val < 0.05) # 1194

##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_gene[,Amyg_Index])$counts, 
	genes = rowData(rse_gene))
dge_Amyg = calcNormFactors(dge_Amyg)
vGene_Amyg = voom(dge_Amyg,mod_Amyg, plot=TRUE)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outGene_Amyg = topTable(eBGene_Amyg,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene_Amyg = outGene_Amyg[rownames(rse_gene),]
sum(outGene_Amyg$adj.P.Val < 0.05) # 2549

g = c("NFKBID", "TNFAIP8L2", "BLNK", "IKBKB")
outGene_Amyg[match(g, outGene_Amyg$Symbol),]
outGene_sACC[match(g, outGene_sACC$Symbol),]
