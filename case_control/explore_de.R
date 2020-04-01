####
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(recount)

## 
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/degradation_rse_BipSeq_BothRegions.rda")
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
				
## add ancestry 
load("../genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## split by brain region
rse_gene_Amygdala = rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sACC = rse_gene[,rse_gene$BrainRegion == "sACC"]
load("../data/degradation_rse_BipSeq_sACC.rda")
load("../data/degradation_rse_BipSeq_Amygdala.rda")
identical(colnames(rse_gene_Amygdala), colnames(cov_rse_amyg)) # TRUE
identical(colnames(rse_gene_sACC), colnames(cov_rse_sacc)) # TRUE

## filter for expression
geneIndex = rowMeans(recount::getRPKM(rse_gene, "Length")) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]
rse_gene_Amygdala = rse_gene_Amygdala[geneIndex,]
rse_gene_sACC = rse_gene_sACC[geneIndex,]

## DEqual objects
load("/dcl01/lieber/ajaffe/lab/degradation_experiments/sACC_RiboZero/sACC_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda")
degradeStats_sACC = degradeStats[rownames(rse_gene),]
load("/dcl01/lieber/ajaffe/lab/degradation_experiments/Amygdala_RiboZero/Amygdala_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda")
degradeStats_Amygdala = degradeStats[rownames(rse_gene),]
load("/dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/rdas/sACC_Plus_Amygdala_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda")
degradeStats_Both = degradeStats[rownames(rse_gene),]
rm(degradeStats)

########################
##  amygdala ###########
boxplot(snpPC1 ~ Dx, data=colData(rse_gene_Amygdala))
boxplot(ERCCsumLogErr ~ Dx, data=colData(rse_gene_Amygdala))
boxplot(rRNA_rate ~ Dx, data=colData(rse_gene_Amygdala))
modAmyg = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene_Amygdala))

##### GENE ######
dgeAmyg = DGEList(counts = assays(rse_gene_Amygdala)$counts, 
	genes = rowData(rse_gene_Amygdala))
dgeAmyg = calcNormFactors(dgeAmyg)
vGeneAmyg = voom(dgeAmyg,modAmyg, plot=TRUE)

fitGeneAmyg = lmFit(vGeneAmyg)
eBGeneAmyg = eBayes(fitGeneAmyg)
outGeneAmyg = topTable(eBGeneAmyg,coef=2,
	p.value = 1,number=nrow(rse_gene_Amygdala))
outGeneAmyg = outGeneAmyg[rownames(rse_gene_Amygdala),]
sum(outGeneAmyg$adj.P.Val < 0.05)

## dequal
plot(outGeneAmyg$t, degradeStats_Amygdala$t)
cor(outGeneAmyg$t, degradeStats_Amygdala$t)

## qSVA - region-specific
degPca_amygOnly = prcomp(t(log2(assays(cov_rse_amyg)$counts+1)))
k_amygOnly = sva::num.sv(log2(assays(cov_rse_amyg)$counts+1), modAmyg) # 16
qSVs_amygOnly = degPca_amygOnly$x[,1:k_amygOnly]

modAmyg_qSVs_amygOnly = cbind(modAmyg,qSVs_amygOnly)
vGeneAmyg_qSVs_amygOnly = voom(dgeAmyg,modAmyg_qSVs_amygOnly, plot=TRUE)
fitGeneAmyg_qSVs_amygOnly = lmFit(vGeneAmyg_qSVs_amygOnly)
eBGeneAmyg_qSVs_amygOnly = eBayes(fitGeneAmyg_qSVs_amygOnly)
outGeneAmyg_qSVs_amygOnly = topTable(eBGeneAmyg_qSVs_amygOnly,coef=2,
	p.value = 1,number=nrow(rse_gene_Amygdala))
outGeneAmyg_qSVs_amygOnly = outGeneAmyg_qSVs_amygOnly[rownames(rse_gene_Amygdala),]
sum(outGeneAmyg_qSVs_amygOnly$adj.P.Val < 0.05)
outGeneAmyg_qSVs_amygOnly$Symbol[outGeneAmyg_qSVs_amygOnly$adj.P.Val < 0.05]

## dequal
plot(outGeneAmyg_qSVs_amygOnly$t, degradeStats_Amygdala$t)
cor(outGeneAmyg_qSVs_amygOnly$t, degradeStats_Amygdala$t)

## qSVA - overall
cov_rse_amygSub = cov_rse[,colnames(rse_gene_Amygdala)]
degPca_amygSub = prcomp(t(log2(assays(cov_rse_amygSub)$counts+1)))
k_amygSub = sva::num.sv(log2(assays(cov_rse_amygSub)$counts+1), modAmyg) # 15
qSVs_amygSub = degPca_amygSub$x[,1:k_amygSub]

modAmyg_qSVs_amygSub = cbind(modAmyg,qSVs_amygSub)
vGeneAmyg_qSVs_amygSub = voom(dgeAmyg,modAmyg_qSVs_amygSub, plot=TRUE)
fitGeneAmyg_qSVs_amygSub = lmFit(vGeneAmyg_qSVs_amygSub)
eBGeneAmyg_qSVs_amygSub = eBayes(fitGeneAmyg_qSVs_amygSub)
outGeneAmyg_qSVs_amygSub = topTable(eBGeneAmyg_qSVs_amygSub,coef=2,
	p.value = 1,number=nrow(rse_gene_Amygdala))
outGeneAmyg_qSVs_amygSub = outGeneAmyg_qSVs_amygSub[rownames(rse_gene_Amygdala),]
sum(outGeneAmyg_qSVs_amygSub$adj.P.Val < 0.05)

plot(outGeneAmyg_qSVs_amygOnly$t, outGeneAmyg_qSVs_amygSub$t)
plot(outGeneAmyg_qSVs_amygOnly$logFC, outGeneAmyg_qSVs_amygSub$logFC)
table(outGeneAmyg_qSVs_amygOnly$adj.P.Val < 0.05, outGeneAmyg_qSVs_amygSub$adj.P.Val < 0.05)

## dequal
plot(outGeneAmyg_qSVs_amygSub$t, degradeStats_Amygdala$t)
cor(outGeneAmyg_qSVs_amygSub$t, degradeStats_Amygdala$t)


########################
##  sACC ###########
modSacc = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene_sACC))

##### GENE ######
dgeSacc = DGEList(counts = assays(rse_gene_sACC)$counts, 
	genes = rowData(rse_gene_sACC))
dgeSacc = calcNormFactors(dgeSacc)
vGeneSacc = voom(dgeSacc,modSacc, plot=TRUE)

fitGeneSacc = lmFit(vGeneSacc)
eBGeneSacc = eBayes(fitGeneSacc)
outGeneSacc = topTable(eBGeneSacc,coef=2,
	p.value = 1,number=nrow(rse_gene_sACC))
outGeneSacc = outGeneSacc[rownames(rse_gene_sACC),]
sum(outGeneSacc$adj.P.Val < 0.05)

plot(outGeneSacc$t, degradeStats_sACC$t)

## qSVA - region-specific
degPca_saccOnly = prcomp(t(log2(assays(cov_rse_sacc)$counts+1)))
k_saccOnly = sva::num.sv(log2(assays(cov_rse_sacc)$counts+1), modSacc) # 14
qSVs_saccOnly = degPca_saccOnly$x[,1:k_saccOnly]

modSacc_qSVs_saccOnly = cbind(modSacc,qSVs_saccOnly)
vGeneSacc_qSVs_saccOnly = voom(dgeSacc,modSacc_qSVs_saccOnly, plot=TRUE)
fitGeneSacc_qSVs_saccOnly = lmFit(vGeneSacc_qSVs_saccOnly)
eBGeneSacc_qSVs_saccOnly = eBayes(fitGeneSacc_qSVs_saccOnly)
outGeneSacc_qSVs_saccOnly = topTable(eBGeneSacc_qSVs_saccOnly,coef=2,
	p.value = 1,number=nrow(rse_gene_sACC))
outGeneSacc_qSVs_saccOnly = outGeneSacc_qSVs_saccOnly[rownames(rse_gene_sACC),]
sum(outGeneSacc_qSVs_saccOnly$adj.P.Val < 0.05)
# outGeneSacc_qSVs_saccOnly$Symbol[outGeneSacc_qSVs_saccOnly$adj.P.Val < 0.05]
## dequal
plot(outGeneSacc_qSVs_saccOnly$t, degradeStats_sACC$t)
cor(outGeneSacc_qSVs_saccOnly$t, degradeStats_sACC$t)

## qSVA - overall
cov_rse_saccSub = cov_rse[,colnames(rse_gene_sACC)]
degPca_saccSub = prcomp(t(log2(assays(cov_rse_saccSub)$counts+1)))
k_saccSub = sva::num.sv(log2(assays(cov_rse_saccSub)$counts+1), modSacc) # 14
qSVs_saccSub = degPca_saccSub$x[,1:k_saccSub]

modSacc_qSVs_saccSub = cbind(modSacc,qSVs_saccSub)
vGeneSacc_qSVs_saccSub = voom(dgeSacc,modSacc_qSVs_saccSub, plot=TRUE)
fitGeneSacc_qSVs_saccSub = lmFit(vGeneSacc_qSVs_saccSub)
eBGeneSacc_qSVs_saccSub = eBayes(fitGeneSacc_qSVs_saccSub)
outGeneSacc_qSVs_saccSub = topTable(eBGeneSacc_qSVs_saccSub,coef=2,
	p.value = 1,number=nrow(rse_gene_sACC))
outGeneSacc_qSVs_saccSub = outGeneSacc_qSVs_saccSub[rownames(rse_gene_sACC),]
sum(outGeneSacc_qSVs_saccSub$adj.P.Val < 0.05)

plot(outGeneSacc_qSVs_saccOnly$t, outGeneSacc_qSVs_saccSub$t)
plot(outGeneSacc_qSVs_saccOnly$logFC, outGeneSacc_qSVs_saccSub$logFC)
table(outGeneSacc_qSVs_saccOnly$adj.P.Val < 0.05, outGeneSacc_qSVs_saccSub$adj.P.Val < 0.05)

table(outGeneSacc_qSVs_saccOnly$adj.P.Val < 0.05, outGeneSacc_qSVs_saccSub$P.Value < 0.01)
table(outGeneSacc_qSVs_saccOnly$P.Value < 0.01, outGeneSacc_qSVs_saccSub$adj.P.Val < 0.05)

#############################
## combined then stratified #
#############################
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

# model
modMain = model.matrix(~Dx + BrainRegion + AgeDeath + Sex + 
	snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

vGeneMain = voom(dge,modMain, plot=TRUE)
fitGeneMain = lmFit(vGeneMain)
eBGeneMain = eBayes(fitGeneMain)
outGeneMain = topTable(eBGeneMain,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneMain = outGeneMain[rownames(rse_gene),]

plot(outGeneMain$t, degradeStats_Both$t)
cor(outGeneMain$t, degradeStats_Both$t)
plot(outGeneMain$t, degradeStats_Both$t_interaction)
cor(outGeneMain$t, degradeStats_Amygdala$t)
cor(outGeneMain$t, degradeStats_sACC$t)

## adding qSVs
degPca = prcomp(t(log2(assays(cov_rse)$counts+1)))
k = sva::num.sv(log2(assays(cov_rse)$counts+1), modMain) # 18
qSVs = degPca$x[,1:k] 

vGeneMain_qSVA = voom(dge,cbind(modMain,qSVs), plot=TRUE)
fitGeneMain_qSVA = lmFit(vGeneMain_qSVA)
eBGeneMain_qSVA = eBayes(fitGeneMain_qSVA)
outGeneMain_qSVA = topTable(eBGeneMain_qSVA,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneMain_qSVA = outGeneMain_qSVA[rownames(rse_gene),]

plot(outGeneMain_qSVA$t, degradeStats_Both$t)
cor(outGeneMain_qSVA$t, degradeStats_Both$t)
plot(outGeneMain_qSVA$t, degradeStats_Both$t_interaction)
cor(outGeneMain_qSVA$t, degradeStats_Amygdala$t)
cor(outGeneMain_qSVA$t, degradeStats_sACC$t)

## interaction
modInt = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + 
	snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

vGeneInt_qSVA = voom(dge,cbind(modInt,qSVs), plot=TRUE)
fitGeneInt_qSVA = lmFit(vGeneInt_qSVA)
eBGeneInt_qSVA = eBayes(fitGeneInt_qSVA)
outGeneInt_qSVA = topTable(eBGeneInt_qSVA,coef=ncol(modInt), # int term
	p.value = 1,number=nrow(rse_gene))
outGeneInt_qSVA = outGeneInt_qSVA[rownames(rse_gene),]
outGeneInt_qSVA_mainTerm = topTable(eBGeneInt_qSVA,coef=2, # dx term
	p.value = 1,number=nrow(rse_gene))
outGeneInt_qSVA_mainTerm = outGeneInt_qSVA_mainTerm[rownames(rse_gene),]

sum(outGeneInt_qSVA$adj.P.Val < 0.2)
sum(outGeneInt_qSVA_mainTerm$adj.P.Val < 0.2)
outGeneInt_qSVA[which(outGeneInt_qSVA$adj.P.Val < 0.2),]
plot(outGeneInt_qSVA$t, degradeStats_Both$t_interaction)

##################################
## posthoc sACC & Amygdala
############################

## sACC
qSVsSub_sacc = qSVs[colnames(rse_gene_sACC),]
modSacc_qSVsSub = cbind(modSacc,qSVsSub_sacc)
vGeneSacc_qSVsSub = voom(dgeSacc,modSacc_qSVsSub, plot=TRUE)
fitGeneSacc_qSVsSub = lmFit(vGeneSacc_qSVsSub)
eBGeneSacc_qSVsSub = eBayes(fitGeneSacc_qSVsSub)
outGeneSacc_qSVsSub = topTable(eBGeneSacc_qSVsSub,coef=2,
	p.value = 1,number=nrow(rse_gene_sACC))
outGeneSacc_qSVsSub = outGeneSacc_qSVsSub[rownames(rse_gene_sACC),]
sum(outGeneSacc_qSVsSub$adj.P.Val < 0.05)
plot(outGeneSacc_qSVsSub$t, outGeneSacc_qSVs_saccOnly$t)
plot(outGeneSacc_qSVsSub$logFC, outGeneSacc_qSVs_saccOnly$logFC)
cor(outGeneSacc_qSVsSub$t, degradeStats_sACC$t)

## Amygdala
qSVsSub_amyg = qSVs[colnames(rse_gene_Amygdala),]
modAmyg_qSVsSub = cbind(modAmyg,qSVsSub_amyg)
vGeneAmyg_qSVsSub = voom(dgeAmyg,modAmyg_qSVsSub, plot=TRUE)
fitGeneAmyg_qSVsSub = lmFit(vGeneAmyg_qSVsSub)
eBGeneAmyg_qSVsSub = eBayes(fitGeneAmyg_qSVsSub)
outGeneAmyg_qSVsSub = topTable(eBGeneAmyg_qSVsSub,coef=2,
	p.value = 1,number=nrow(rse_gene_sACC))
outGeneAmyg_qSVsSub = outGeneAmyg_qSVsSub[rownames(rse_gene_sACC),]
sum(outGeneAmyg_qSVsSub$adj.P.Val < 0.05)
plot(outGeneAmyg_qSVsSub$t, outGeneAmyg_qSVs_amygOnly$t)
plot(outGeneAmyg_qSVsSub$logFC, outGeneAmyg_qSVs_amygOnly$logFC)
cor(outGeneAmyg_qSVsSub$t, degradeStats_Amygdala$t)

## cross region compare
plot(outGeneSacc$t, outGeneAmyg$t)
plot(outGeneSacc_qSVsSub$t, outGeneAmyg_qSVsSub$t)
plot(outGeneSacc_qSVs_saccSub$t, outGeneAmyg_qSVs_amygSub$t)

table(outGeneSacc$adj.P.Val < 0.05, outGeneAmyg$adj.P.Val < 0.05)
table(outGeneSacc_qSVsSub$adj.P.Val < 0.05, outGeneAmyg_qSVsSub$adj.P.Val < 0.05)
table(outGeneSacc_qSVs_saccSub$adj.P.Val < 0.05, outGeneAmyg_qSVs_amygSub$adj.P.Val < 0.05)