##### 
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)
library(purrr)
library(here)

## load data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
# load("../data/zandiHypde_bipolar_rseExon_n511.rda")
# load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
# load("../data/zandiHypde_bipolar_rseTx_n511.rda")

load(here("data","degradation_rse_BipSeq_BothRegions.rda"), verbose = TRUE)

## deconvolution results
load(here("deconvolution","est_prop_Bisque.Rdata"),verbose = TRUE)
# # add counts to Tx rse
# load("../data/zandiHypde_bipolar_txCounts_n511.rda")
# identical(colnames(txNumReads), colData(rse_tx)$SAMPLE_ID)  ## TRUE
# colnames(txNumReads) = rownames(colData(rse_tx))
# assays(rse_tx)$counts = txNumReads

identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
identical(colnames(rse_gene), rownames(est_prop_bisque$bulk.props)) #TRUE

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

# ## exon
# assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
# exonIndex = rowMeans(assays(rse_exon)$rpkm) > 0.3
# rse_exon = rse_exon[exonIndex,]
# 
# ## junction
# rowRanges(rse_jxn)$Length <- 100
# assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
# jxnIndex = rowMeans(assays(rse_jxn)$rp10m) > 0.35 & rowData(rse_jxn)$Class != "Novel"
# rse_jxn = rse_jxn[jxnIndex,]
# 
# ## transcript
# txIndex = rowMeans(assays(rse_tx)$tpm) > 0.4 
# rse_tx = rse_tx[txIndex,]

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
modSep = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

modSep_prop <- cbind(modSep,est_prop_bisque$bulk.props[,c("Astro","Micro","Oligo","OPC","Excit")])
modSep_ilr <- cbind(modSep,est_prop_bisque$ilr)

modSep_deconvo <- list(prop = modSep_prop, ilr = modSep_ilr)
#########################
## split back by region #
#########################

#### both ####
sACC_Index = which(colData(rse_gene)$BrainRegion == "sACC")
mod_sACC <- map(modSep_deconvo, ~cbind(.x[sACC_Index,], qSV_mat[sACC_Index, ]))

Amyg_Index = which(colData(rse_gene)$BrainRegion == "Amygdala")
mod_Amyg <- map(modSep_deconvo, ~cbind(.x[Amyg_Index,], qSV_mat[Amyg_Index, ]))

#################
### Gene ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_gene[,sACC_Index])$counts, 
                   genes = rowData(rse_gene))
dge_sACC = calcNormFactors(dge_sACC)

outGene_sACC <- map(mod_sACC, function(m){

  vGene_sACC = voom(dge_sACC,m, plot=FALSE)
  fitGene_sACC = lmFit(vGene_sACC)
  eBGene_sACC = eBayes(fitGene_sACC)
  outGene_sACC = topTable(eBGene_sACC,coef=2,
                          p.value = 1,number=nrow(rse_gene))
  outGene_sACC = outGene_sACC[rownames(rse_gene),]
  return(outGene_sACC)
})

map_int(outGene_sACC, ~sum(.x$adj.P.Val < 0.05))
# prop  ilr 
# 318  379 

##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_gene[,Amyg_Index])$counts, 
	genes = rowData(rse_gene))
dge_Amyg = calcNormFactors(dge_Amyg)

outGene_Amyg <- map(mod_Amyg, function(m){

  vGene_Amyg = voom(dge_Amyg,m, plot=FALSE)
  
  fitGene_Amyg = lmFit(vGene_Amyg)
  eBGene_Amyg = eBayes(fitGene_Amyg)
  outGene_Amyg = topTable(eBGene_Amyg,coef=2,
                          p.value = 1,number=nrow(rse_gene))
  outGene_Amyg = outGene_Amyg[rownames(rse_gene),]
  return(outGene_Amyg)
})

map_int(outGene_Amyg, ~sum(.x$adj.P.Val < 0.05))
# prop  ilr 
# 324  258 
# #################
# ### Exon ########
# #################
# 
# ##### sACC ######
# dee_sACC = DGEList(counts = assays(rse_exon[,sACC_Index])$counts, 
# 	genes = rowData(rse_exon))
# dee_sACC = calcNormFactors(dee_sACC)
# vExon_sACC = voom(dee_sACC,mod_sACC, plot=TRUE)
# 
# fitExon_sACC = lmFit(vExon_sACC)
# eBExon_sACC = eBayes(fitExon_sACC)
# outExon_sACC = topTable(eBExon_sACC,coef=2,
# 	p.value = 1,number=nrow(rse_exon))
# outExon_sACC = outExon_sACC[rownames(rse_exon),]
# sum(outExon_sACC$adj.P.Val < 0.05)
# 
# ##### Amygdala ######
# dee_Amyg = DGEList(counts = assays(rse_exon[,Amyg_Index])$counts, 
# 	genes = rowData(rse_exon))
# dee_Amyg = calcNormFactors(dee_Amyg)
# vExon_Amyg = voom(dee_Amyg,mod_Amyg, plot=TRUE)
# 
# fitExon_Amyg = lmFit(vExon_Amyg)
# eBExon_Amyg = eBayes(fitExon_Amyg)
# outExon_Amyg = topTable(eBExon_Amyg,coef=2,
# 	p.value = 1,number=nrow(rse_exon))
# outExon_Amyg = outExon_Amyg[rownames(rse_exon),]
# sum(outExon_Amyg$adj.P.Val < 0.05)
# 
# 
# #################
# ### Junction ########
# #################
# 
# ##### sACC ######
# dje_sACC = DGEList(counts = assays(rse_jxn[,sACC_Index])$counts, 
# 	genes = rowData(rse_jxn))
# dje_sACC = calcNormFactors(dje_sACC)
# vJxn_sACC = voom(dje_sACC,mod_sACC, plot=TRUE)
# 
# fitJxn_sACC = lmFit(vJxn_sACC)
# eBJxn_sACC = eBayes(fitJxn_sACC)
# outJxn_sACC = topTable(eBJxn_sACC,coef=2,
# 	p.value = 1,number=nrow(rse_jxn))
# outJxn_sACC = outJxn_sACC[rownames(rse_jxn),]
# sum(outJxn_sACC$adj.P.Val < 0.05)
# 
# ##### Amygdala ######
# dje_Amyg = DGEList(counts = assays(rse_jxn[,Amyg_Index])$counts, 
# 	genes = rowData(rse_jxn))
# dje_Amyg = calcNormFactors(dje_Amyg)
# vJxn_Amyg = voom(dje_Amyg,mod_Amyg, plot=TRUE)
# 
# fitJxn_Amyg = lmFit(vJxn_Amyg)
# eBJxn_Amyg = eBayes(fitJxn_Amyg)
# outJxn_Amyg = topTable(eBJxn_Amyg,coef=2,
# 	p.value = 1,number=nrow(rse_jxn))
# outJxn_Amyg = outJxn_Amyg[rownames(rse_jxn),]
# sum(outJxn_Amyg$adj.P.Val < 0.05)
# 
# 
# #################
# ### Transcript ########
# #################
# 
# txExprs = log2(assays(rse_tx)$tpm+ 1)
# 
# ##### sACC ######
# # fitTx_sACC = lmFit(txExprs[,sACC_Index], modSep[sACC_Index,])
# fitTx_sACC = lmFit(txExprs[,sACC_Index], mod_sACC)
# eBTx_sACC = eBayes(fitTx_sACC)
# outTx_sACC = topTable(eBTx_sACC,coef=2,
# 	p.value = 1,number=nrow(rse_tx), 
# 	genelist = rowRanges(rse_tx))
# outTx_sACC = outTx_sACC[rownames(rse_tx),c(28:33, 10:11, 13, 15:17, 19, 26)]
# sum(outTx_sACC$adj.P.Val < 0.05)  ## 99
# 
# ##### Amygdala ######
# # fitTx_Amyg = lmFit(txExprs[,Amyg_Index], modSep[Amyg_Index,])
# fitTx_Amyg = lmFit(txExprs[,Amyg_Index], mod_Amyg)
# eBTx_Amyg = eBayes(fitTx_Amyg)
# outTx_Amyg = topTable(eBTx_Amyg,coef=2,
# 	p.value = 1,number=nrow(rse_tx),
# 	genelist = rowRanges(rse_tx))
# outTx_Amyg = outTx_Amyg[rownames(rse_tx),c(28:33, 10:11, 13, 15:17, 19, 26)]
# sum(outTx_Amyg$adj.P.Val < 0.05)  ## 3
# 


####################
### core output ####
nam = c("logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B", "deconvo_terms")

outGene_sACC_merged <- do.call("rbind",outGene_sACC)
outGene_sACC_merged$deconvo_terms <- ss(rownames(outGene_sACC_merged),"\\.")

outGene_Amyg_merged <- do.call("rbind",outGene_Amyg)
outGene_Amyg_merged$deconvo_terms <- ss(rownames(outGene_Amyg_merged),"\\.")

geneOut = cbind(outGene_Amyg_merged[,nam], outGene_sACC_merged[,nam])
colnames(geneOut) = paste0(colnames(geneOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
geneOut_deconvo <- geneOut

# exonOut = cbind(outExon_Amyg[,nam], outExon_sACC[,nam])
# colnames(exonOut) = paste0(colnames(exonOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
# 
# jxnOut = cbind(outJxn_Amyg[,nam], outJxn_sACC[,nam])
# colnames(jxnOut) = paste0(colnames(jxnOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
# 
# txOut = cbind(outTx_Amyg[,nam], outTx_sACC[,nam])
# colnames(txOut) = paste0(colnames(txOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
# 
# statOut = rbind(geneOut, exonOut, jxnOut, txOut)
# save(statOut, compress=TRUE, file = "bipolarControl_deStats_byRegion_qSVAjoint.rda")
save(geneOut_deconvo, file = "bipolarControl_deStats_byRegion_qSVAjoint_deconvo.rda")
# ######################################
# ##### interaction/cross-region #######

## overall model
mod = cbind(modJoint, qSV_mat)

mod_prop <- cbind(mod ,est_prop_bisque$bulk.props[,c("Astro","Micro","Oligo","OPC","Excit")])
mod_ilr <- cbind(mod ,est_prop_bisque$ilr)

mod_deconvo <- list(prop = modSep_prop, ilr = modSep_ilr)
# ###########
# ## Gene ###
# ###########

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

eBGene <- map(mod_deconvo,function(mod){
  vGene = voom(dge,mod, plot=FALSE)
  
  ## do duplicate correlation
  gene_dupCorr = duplicateCorrelation(vGene$E, mod, block=colData(rse_gene)$BrNum)
  save(gene_dupCorr, file = "geneLevel_duplicateCorrelation_limma_forDE.rda")
  load("geneLevel_duplicateCorrelation_limma_forDE.rda")
  
  # and then fit
  fitGene = lmFit(vGene, mod,
                  correlation=gene_dupCorr$consensus.correlation,
                  block=colData(rse_gene)$BrNum)
  eBGene = eBayes(fitGene)
  return(eBGene)
})

outGene_mainEffect = map(eBGene, ~topTable(.x,coef=2, p.value = 1,number=nrow(rse_gene), sort="none"))

outGene_interactionEffect = map(eBGene, ~topTable(.x, coef=ncol(modJoint),p.value = 1,number=nrow(rse_gene), sort="none"))

map_int(outGene_mainEffect, ~sum(.x$adj.P.Val < 0.05))
map_int(outGene_interactionEffect, ~sum(.x$adj.P.Val < 0.05))

# ###########
# ## Exon ###
# ###########
# 
# dee = DGEList(counts = assays(rse_exon)$counts, 
# 	genes = rowData(rse_exon))
# dee = calcNormFactors(dee)
# vExon = voom(dee,mod, plot=TRUE)
# 
# ## do duplicate correlation
# exon_dupCorr = duplicateCorrelation(vExon$E, mod, block=colData(rse_exon)$BrNum)
# save(exon_dupCorr, file = "exonLevel_duplicateCorrelation_limma_forDE.rda")
# load("exonLevel_duplicateCorrelation_limma_forDE.rda")
# 
# # and then fit
# fitExon = lmFit(vExon, mod, 
# 	correlation=exon_dupCorr$consensus.correlation, 
# 	block=colData(rse_exon)$BrNum)
# eBExon = eBayes(fitExon)
# outExon_mainEffect = topTable(eBExon,coef=2,
# 	p.value = 1,number=nrow(rse_exon), sort="none")
# outExon_interactionEffect = topTable(eBExon,coef=ncol(modJoint),
# 	p.value = 1,number=nrow(rse_exon), sort="none")
# 
# sum(outExon_mainEffect$adj.P.Val < 0.05)
# sum(outExon_interactionEffect$adj.P.Val < 0.05)
# 
# 
# ###########
# ## Junction ###
# ###########
# 
# dje = DGEList(counts = assays(rse_jxn)$counts, 
# 	genes = rowData(rse_jxn))
# dje = calcNormFactors(dje)
# vJxn = voom(dje,mod, plot=TRUE)
# 
# ## do duplicate correlation
# # jxn_dupCorr = duplicateCorrelation(vJxn$E, mod, block=colData(rse_jxn)$BrNum)
# # save(jxn_dupCorr, file = "jxnLevel_duplicateCorrelation_limma_forDE.rda")
# load("jxnLevel_duplicateCorrelation_limma_forDE.rda")
# # and then fit
# fitJxn = lmFit(vJxn, mod,
# 	correlation=jxn_dupCorr$consensus.correlation, 
# 	block=colData(rse_jxn)$BrNum)
# eBJxn = eBayes(fitJxn)
# outJxn_mainEffect = topTable(eBJxn,coef=2,
# 	p.value = 1,number=nrow(rse_jxn), sort="none")
# outJxn_interactionEffect = topTable(eBJxn,coef=ncol(modJoint),
# 	p.value = 1,number=nrow(rse_jxn), sort="none")
# 
# sum(outJxn_mainEffect$adj.P.Val < 0.05)
# sum(outJxn_interactionEffect$adj.P.Val < 0.05)
# 
# ###########
# ## Txs ###
# ###########
# 
# ## do duplicate correlation
# # tx_dupCorr = duplicateCorrelation(txExprs, mod, block=colData(rse_tx)$BrNum)
# # save(tx_dupCorr, file = "txLevel_duplicateCorrelation_limma_forDE.rda")
# load("txLevel_duplicateCorrelation_limma_forDE.rda")
# 
# # and then fit
# fitTx = lmFit(txExprs, mod, 
# 	correlation=tx_dupCorr$consensus.correlation, 
# 	block=colData(rse_tx)$BrNum)
# eBTx = eBayes(fitTx)
# outTx_mainEffect = topTable(eBTx,coef=2,
# 	p.value = 1,number=nrow(rse_tx), sort="none")
# outTx_interactionEffect = topTable(eBTx,coef=ncol(modJoint),
# 	p.value = 1,number=nrow(rse_tx), sort="none")
# 
# sum(outTx_mainEffect$adj.P.Val < 0.05)
# sum(outTx_interactionEffect$adj.P.Val < 0.05)
# 
# ######################
# ## save outputs ######
# ######################
# 
# outGene_bothRegion = outGene_mainEffect
# colnames(outGene_bothRegion)[c(11,13:16)] = paste0(colnames(outGene_bothRegion)[c(11,13:16)],"_dxEffect")
# outGene_bothRegion = cbind(outGene_bothRegion, outGene_interactionEffect[,c(11,13:16)])
# colnames(outGene_bothRegion)[17:21] = paste0(colnames(outGene_bothRegion)[17:21],"_intEffect")
# outGene_bothRegion = outGene_bothRegion[,c(5,2,11,13:21, 3:4,6:10,12,1)]
# 
# outExon_bothRegion = outExon_mainEffect
# colnames(outExon_bothRegion)[c(11,13:16)] = paste0(colnames(outExon_bothRegion)[c(11,13:16)],"_dxEffect")
# outExon_bothRegion = cbind(outExon_bothRegion, outExon_interactionEffect[,c(11,13:16)])
# colnames(outExon_bothRegion)[17:21] = paste0(colnames(outExon_bothRegion)[17:21],"_intEffect")
# outExon_bothRegion = outExon_bothRegion[,c(5,2,11,13:21, 3:4,6:10,12,1)]
# 
# outJxn_bothRegion = outJxn_mainEffect
# colnames(outJxn_bothRegion)[c(18,20:23)] = paste0(colnames(outJxn_bothRegion)[c(18,20:23)],"_dxEffect")
# outJxn_bothRegion = cbind(outJxn_bothRegion, outJxn_interactionEffect[,c(18,20:23)])
# colnames(outJxn_bothRegion)[24:28] = paste0(colnames(outJxn_bothRegion)[24:28],"_intEffect")
# outJxn_bothRegion = outJxn_bothRegion[,c(14,13, 18,20:28, 1:12, 15, 19, 17)]
# 
# outTx_bothRegion = cbind(rowData(rse_tx), outTx_mainEffect)
# colnames(outTx_bothRegion)[c(23,25:28)] = paste0(colnames(outTx_bothRegion)[c(23,25:28)],"_dxEffect")
# outTx_bothRegion = cbind(outTx_bothRegion, outTx_interactionEffect[,-2])
# colnames(outTx_bothRegion)[29:33] = paste0(colnames(outTx_bothRegion)[29:33],"_intEffect")
# outTx_bothRegion = outTx_bothRegion[,c(8,5,23,25:33, 24,1:4, 6:7,9:22)]
# 
# save(outGene_bothRegion, outExon_bothRegion, outJxn_bothRegion, outTx_bothRegion,
# 	file = "interaction_model_results.rda", compress=TRUE)
# 
# ### write csv
# dir.create("csv")
# outGene_bothRegion$gencodeTx = NULL
# write.csv(outGene_bothRegion, file = gzfile("csv/geneLevel_jointModeling_lmer.csv.gz"))
# outExon_bothRegion$gencodeTx = NULL
# write.csv(outExon_bothRegion, file = gzfile("csv/exonLevel_jointModeling_lmer.csv.gz"))
# outJxn_bothRegion$gencodeTx = NULL
# write.csv(outJxn_bothRegion, file = gzfile("csv/jxnLevel_jointModeling_lmer.csv.gz"))
# outTx_bothRegion$gencodeTx = NULL
# write.csv(outTx_bothRegion, file = gzfile("csv/txLevel_jointModeling_lmer.csv.gz"))
# 
# ## plots
# smoothScatter(statOut$logFC_Amyg, statOut$logFC_sACC,pch = 21, bg="grey",
# 	xlim = c(-1,1), ylim = c(-1,1))
# gIndex = grep("ENSG", rownames(statOut))
# smoothScatter(statOut$logFC_Amyg[gIndex],
# 	statOut$logFC_sACC[gIndex],pch = 21, bg="grey",
# 	xlim = c(-1,1), ylim = c(-1,1))
# plot(statOut$logFC_Amyg[gIndex], statOut$logFC_sACC[gIndex],
# 	pch = 21, bg="grey",xlim = c(-1,1), ylim = c(-1,1),cex=0.2)
# 	
# jIndex = grep("chr", rownames(statOut))
# smoothScatter(statOut$logFC_Amyg[jIndex],
# 	statOut$logFC_sACC[jIndex],xlim = c(-1,1), ylim = c(-1,1))
