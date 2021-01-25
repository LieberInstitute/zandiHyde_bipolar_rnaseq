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

###################
## add drug info ##
###################

drug_info = read.delim("../drug_info.tsv",as.is=TRUE)
rse_gene$lithium = drug_info$lithium[match(rse_gene$BrNum, drug_info$BrNum)]
rse_gene$lifetime_lithium = drug_info$lifetime_lithium[match(rse_gene$BrNum, drug_info$BrNum)]
rse_gene$lithium_group = NA
rse_gene$lithium_group[which(rse_gene$lithium)] =1 
rse_gene$lithium_group[which(!rse_gene$lithium & !rse_gene$lifetime_lithium)] =0

# model w/o interaction to subset by region
modSep_lithium = model.matrix(~lithium_group + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))
rse_gene_lithium = rse_gene[,rownames(modSep_lithium)]
qSV_mat_lithium = qSV_mat[rownames(modSep_lithium),]

#########################
## split back by region #
#########################

#### both ####
sACC_Index = which(colData(rse_gene_lithium)$BrainRegion == "sACC")
mod_sACC = cbind(modSep_lithium[sACC_Index,], qSV_mat_lithium[sACC_Index, ])
Amyg_Index = which(colData(rse_gene_lithium)$BrainRegion == "Amygdala")
mod_Amyg = cbind(modSep_lithium[Amyg_Index,], qSV_mat_lithium[Amyg_Index, ])

#################
### Gene ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_gene_lithium[,sACC_Index])$counts, 
	genes = rowData(rse_gene_lithium))
dge_sACC = calcNormFactors(dge_sACC)
vGene_sACC = voom(dge_sACC,mod_sACC, plot=TRUE)

fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outGene_sACC = topTable(eBGene_sACC,coef=2,
	p.value = 1,number=nrow(rse_gene_lithium))
outGene_sACC = outGene_sACC[rownames(rse_gene_lithium),]
sum(outGene_sACC$adj.P.Val < 0.05)
min(outGene_sACC$adj.P.Val)

##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_gene_lithium[,Amyg_Index])$counts, 
	genes = rowData(rse_gene_lithium))
dge_Amyg = calcNormFactors(dge_Amyg)
vGene_Amyg = voom(dge_Amyg,mod_Amyg, plot=TRUE)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outGene_Amyg = topTable(eBGene_Amyg,coef=2,
	p.value = 1,number=nrow(rse_gene_lithium))
outGene_Amyg = outGene_Amyg[rownames(rse_gene_lithium),]
sum(outGene_Amyg$adj.P.Val < 0.05)

### core output ####
nam = c("logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")
geneOut = cbind(outGene_Amyg[,nam], outGene_sACC[,nam])
colnames(geneOut) = paste0(colnames(geneOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
geneOut$Symbol = rowData(rse_gene)$Symbol
geneOut = geneOut[,c(13,1:12)]
write.csv(geneOut, file = "lithium_effects_bothRegions.csv")
###################################
## read back in bipolar effects ###

load("bipolarControl_deStats_byRegion_qSVAjoint.rda",verbose=TRUE)
bipolar_stats= statOut[rownames(rse_gene_lithium),]

## Amygdala ##
## all genes 
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main = 2)
plot(bipolar_stats$t_Amyg, outGene_Amyg$t,
	pch=21,bg="grey", main = "Amygdala (T-stats)",
	xlab=  "BPD vs Control Effect", ylab = "Lithium Effect")

## bipolar significant
amyg_index =which(bipolar_stats$adj.P.Val_Amyg < 0.05)
plot(bipolar_stats$t_Amyg[amyg_index], outGene_Amyg$t[amyg_index],
	pch=21,bg="grey", main = "Amygdala (T-stats)",
	xlab=  "BPD vs Control Effect", ylab = "Lithium Effect")

	
### sacc
plot(bipolar_stats$t_sACC , outGene_sACC$t,	
	pch=21,bg="grey", main = "sACC (T-stats)",
	xlab=  "BPD vs Control Effect", ylab = "Lithium Effect")
sacc_index =which(bipolar_stats$adj.P.Val_sACC< 0.05)
plot(bipolar_stats$t_sACC[sacc_index] , outGene_sACC$t[sacc_index],	
	pch=21,bg="grey", main = "sACC (T-stats)",
	xlab=  "BPD vs Control Effect", ylab = "Lithium Effect")
abline(h=0,lty=2)