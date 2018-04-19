##### 
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)

## load data
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../data/degradation_rse_BipSeq_BothRegions.rda")

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
## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

## exon
assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
exonIndex = rowMeans(assays(rse_exon)$rpkm) > 0.3
rse_exon = rse_exon[exonIndex,]

## junction
rowRanges(rse_jxn)$Length <- 100
assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
jxnIndex = rowMeans(assays(rse_jxn)$rp10m) > 0.35 & rowData(rse_jxn)$Class != "Novel"
rse_jxn = rse_jxn[jxnIndex,]

## transcript
txIndex = rowMeans(assays(rse_tx)$tpm) > 0.4 
rse_tx = rse_tx[txIndex,]

##############
## get qSVs ##
##############

modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 18
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 87%

# model w/o interaction to subset by region
modSep = modJoint = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

#########################
## split back by region #
#########################

#### both ####
sACC_Index = which(colData(rse_gene)$BrainRegion == "sACC")
mod_sACC = cbind(modSep[sACC_Index,], qSV_mat[sACC_Index, ])
Amyg_Index = which(colData(rse_gene)$BrainRegion == "Amygdala")
mod_Amyg = cbind(modSep[Amyg_Index,], qSV_mat[Amyg_Index, ])

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
sum(outGene_sACC$adj.P.Val < 0.05)

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
sum(outGene_Amyg$adj.P.Val < 0.05)


#################
### Exon ########
#################

##### sACC ######
dee_sACC = DGEList(counts = assays(rse_exon[,sACC_Index])$counts, 
	genes = rowData(rse_exon))
dee_sACC = calcNormFactors(dee_sACC)
vExon_sACC = voom(dee_sACC,mod_sACC, plot=TRUE)

fitExon_sACC = lmFit(vExon_sACC)
eBExon_sACC = eBayes(fitExon_sACC)
outExon_sACC = topTable(eBExon_sACC,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon_sACC = outExon_sACC[rownames(rse_exon),]
sum(outExon_sACC$adj.P.Val < 0.05)

##### Amygdala ######
dee_Amyg = DGEList(counts = assays(rse_exon[,Amyg_Index])$counts, 
	genes = rowData(rse_exon))
dee_Amyg = calcNormFactors(dee_Amyg)
vExon_Amyg = voom(dee_Amyg,mod_Amyg, plot=TRUE)

fitExon_Amyg = lmFit(vExon_Amyg)
eBExon_Amyg = eBayes(fitExon_Amyg)
outExon_Amyg = topTable(eBExon_Amyg,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon_Amyg = outExon_Amyg[rownames(rse_exon),]
sum(outExon_Amyg$adj.P.Val < 0.05)


#################
### Junction ########
#################

##### sACC ######
dje_sACC = DGEList(counts = assays(rse_jxn[,sACC_Index])$counts, 
	genes = rowData(rse_jxn))
dje_sACC = calcNormFactors(dje_sACC)
vJxn_sACC = voom(dje_sACC,mod_sACC, plot=TRUE)

fitJxn_sACC = lmFit(vJxn_sACC)
eBJxn_sACC = eBayes(fitJxn_sACC)
outJxn_sACC = topTable(eBJxn_sACC,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn_sACC = outJxn_sACC[rownames(rse_jxn),]
sum(outJxn_sACC$adj.P.Val < 0.05)

##### Amygdala ######
dje_Amyg = DGEList(counts = assays(rse_jxn[,Amyg_Index])$counts, 
	genes = rowData(rse_jxn))
dje_Amyg = calcNormFactors(dje_Amyg)
vJxn_Amyg = voom(dje_Amyg,mod_Amyg, plot=TRUE)

fitJxn_Amyg = lmFit(vJxn_Amyg)
eBJxn_Amyg = eBayes(fitJxn_Amyg)
outJxn_Amyg = topTable(eBJxn_Amyg,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn_Amyg = outJxn_Amyg[rownames(rse_jxn),]
sum(outJxn_Amyg$adj.P.Val < 0.05)


#################
### Transcript ########
#################

txExprs = log2(assays(rse_tx)$tpm+ 1)

##### sACC ######
fitTx_sACC = lmFit(txExprs[,sACC_Index], modSep[sACC_Index,])
eBTx_sACC = eBayes(fitTx_sACC)
outTx_sACC = topTable(eBTx_sACC,coef=2,
	p.value = 1,number=nrow(rse_tx), 
	genelist = rowRanges(rse_tx))
outTx_sACC = outTx_sACC[rownames(rse_tx),c(28:33, 10:11, 13, 15:17, 19, 26)]
sum(outTx_sACC$adj.P.Val < 0.05)

##### Amygdala ######
fitTx_Amyg = lmFit(txExprs[,Amyg_Index], modSep[Amyg_Index,])
eBTx_Amyg = eBayes(fitTx_Amyg)
outTx_Amyg = topTable(eBTx_Amyg,coef=2,
	p.value = 1,number=nrow(rse_tx),
	genelist = rowRanges(rse_tx))
outTx_Amyg = outTx_Amyg[rownames(rse_tx),c(28:33, 10:11, 13, 15:17, 19, 26)]
sum(outTx_Amyg$adj.P.Val < 0.05)

####################
### core output ####
nam = c("logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")

geneOut = cbind(outGene_Amyg[,nam], outGene_sACC[,nam])
colnames(geneOut) = paste0(colnames(geneOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))

exonOut = cbind(outExon_Amyg[,nam], outExon_sACC[,nam])
colnames(exonOut) = paste0(colnames(exonOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))

jxnOut = cbind(outJxn_Amyg[,nam], outJxn_sACC[,nam])
colnames(jxnOut) = paste0(colnames(jxnOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))

txOut = cbind(outTx_Amyg[,nam], outTx_sACC[,nam])
colnames(txOut) = paste0(colnames(txOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))

statOut = rbind(geneOut, exonOut, jxnOut, txOut)
save(statOut, compress=TRUE, file = "bipolarControl_deStats_byRegion_qSVAjoint.rda")

######################################
##### interaction/cross-region #######

## overall model
mod = cbind(modJoint, qSV_mat)

###########
## Gene ###
###########

dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)
vGene = voom(dge,mod, plot=TRUE)

## do duplicate correlation
gene_dupCorr = duplicateCorrelation(vGene$E, mod, block=colData(rse_gene)$BrNum)
save(gene_dupCorr, file = "geneLevel_duplicateCorrelation_limma_forDE.rda")

# and then fit
fitGene = lmFit(vGene,block = colData(rse_gene)$BrNum, correlation=gene_dupCorr)
eBGene = eBayes(fitGene)
outGene_mainEffect = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene_mainEffect = outGene_mainEffect[rownames(rse_gene),]

outGene_interactionEffect = topTable(eBGene,coef=ncol(mod),
	p.value = 1,number=nrow(rse_gene), sort="none")
