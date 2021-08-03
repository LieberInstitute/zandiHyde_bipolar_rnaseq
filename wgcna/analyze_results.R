##
library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
library(clusterProfiler)
library(readxl)
library(RColorBrewer)

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

##########
## model #
##########
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

## qsva
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 18
qSV_mat = prcomp(t(degExprs))$x[,1:k]

## join and move around region, dx and interaction for cleaning
modQsva = cbind(modJoint[,c(1:4,14,5:13)], qSV_mat)

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

######################
## gene ontology
gList = split(rowData(rse_gene)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene)$EntrezID
univ = as.character(univ[!is.na(univ)])

go = compareCluster(gList, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
save(go, file = "go_enrichment_bpdControl_wgcna.rda")

goDf = as.data.frame(go)
goCheck = goDf[goDf$Cluster %in% c("red", "pink", "magenta", "royalblue") &
	goDf$qvalue < 0.05,]
goCheck  =goCheck[order(goCheck$pvalue),]
write.csv(goCheck, file = "go_enrichment_bpdControl_wgcna_fourCandidateModules.csv")
##############
## associate eigengenes with brain region
m = modQsva[,1:5] # this is what was protected

MEs = net$MEs
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col]

## check
statList = lapply(MEs, function(x) summary(lmer(x ~ m + (1|rse_gene$BrNum) - 1))$coef)

# modified
bpdEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
regionEffect = as.data.frame(t(sapply(statList, function(x) x[3,])))
ageEffect = as.data.frame(t(sapply(statList, function(x) x[4,])))
intEffect = as.data.frame(t(sapply(statList, function(x) x[5,])))
colnames(bpdEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
	"slope", "se", "df", "t", "pvalue")

signif(bpdEffect, 3)
signif(regionEffect, 3)
signif(ageEffect, 3)
signif(intEffect, 3)

##################
# make boxplots ##
lab = paste0(substr(rse_gene$BrainRegion,1,4), ":", rse_gene$Dx)
lab = factor(lab, levels = c("sACC:Control", "sACC:Bipolar", "Amyg:Control", "Amyg:Bipolar"))

pdf("MEs_vs_dx.pdf",w=8,h=6)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=2,cex.lab=1.8,cex.main = 1.8)
for(i in 1:ncol(MEs)) {
	boxplot(MEs[,i] ~ lab, outline = FALSE, xlab="",
		ylim = quantile(unlist(MEs),c(0.002,0.998)),main = colnames(MEs)[i],
		names = gsub(":", "\n", levels(lab)), ylab = "Module Eigengene")
	points(MEs[,i] ~ jitter(as.numeric(lab),amount=0.1), pch=21, bg=lab)
	legend("top", c(paste0("Region p=", signif(regionEffect[i,5],3)), 
		paste0("Dx p=", signif(bpdEffect[i,5],3))),cex=1.4)
}
dev.off()
	
clean_ME = t(cleaningY(t(MEs), m, P=2))
pdf("clean_MEs_vs_dx.pdf",w=5,h=6)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=2,cex.lab=1.8,cex.main = 1.8)
for(i in 1:ncol(clean_ME)) {
	boxplot(clean_ME[,i] ~ rse_gene$Dx, outline = FALSE, xlab="",
		ylim = quantile(unlist(clean_ME),c(0.002,0.998)),main = colnames(MEs)[i],
		ylab = "Module Eigengene (Adj)")
	points(clean_ME[,i] ~ jitter(as.numeric(rse_gene$Dx),amount=0.1), pch=21, bg=lab)
	legend("top", paste0("Dx p=", signif(bpdEffect[i,5],3)),cex=1.4)
}
dev.off()
	

##############################
## enrichment of DEGs #####
#############################
load("../case_control/interaction_model_results.rda")
rm(outExon_bothRegion, outJxn_bothRegion, outTx_bothRegion)

identical(rownames(outGene_bothRegion), fNames) # TRUE
tt = table(net$colorsLab, outGene_bothRegion$adj.P.Val_dxEffect < 0.05)
tt = tt[rownames(bpdEffect),]

plot(prop.table(tt,1)[,2], -log10(bpdEffect$pvalue))

# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
	tab = table(net$colorsLab == cc,outGene_bothRegion$adj.P.Val_dxEffect < 0.05)
	c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

plot(-log10(deModEnrich$Pvalue), -log10(bpdEffect$pvalue ))

## do stuff on red

## read in magma
magma = read_excel("../PGC_BP_Magma_table.xlsx", skip = 2)
magma = as.data.frame(magma)
magma$FDR_JOINT = p.adjust(magma$P_JOINT, "fdr")
magma$BONF_JOINT = p.adjust(magma$P_JOINT, "bonf")
magma$Module = net$colorsLab[match(magma$GENE, rowData(rse_gene)$EntrezID)]
magma$isSig = factor(ifelse(magma$BONF_JOINT < 0.05, "Yes","No"))
ms = unique(magma$Module)
ms = ms[!is.na(ms)]

## test
magmaEnrich = t(sapply(ms, function(m) {
	tt =table(magma$Module == m, magma$isSig)
	c(getOR(tt), chisq.test(tt)$p.value)
}))
colnames(magmaEnrich) = c("OR", "pvalue")
magmaEnrich = as.data.frame(magmaEnrich)

## add counts
magmaEnrich$MagmaSig = table(magma$Module[magma$BONF_JOINT < 0.05])[rownames(magmaEnrich)]
magmaEnrich$NumGenes = table(magma$Module)[rownames(magmaEnrich)]
write.csv(magmaEnrich, file = "magma_module_enrichment.csv")

magmaSig = magma[magma$BONF_JOINT < 0.05 & magma$GENE %in% rowData(rse_gene)$EntrezID,]

entrezIDs = split(rowData(rse_gene)$EntrezID, net$colorsLab)
entrezOverlap = sapply(entrezIDs, function(x) sum(x %in% magmaSig$GENE,na.rm=TRUE))
entrezPresent = sapply(entrezIDs, function(x) sum(x %in% magma$GENE,na.rm=TRUE))
x = data.frame(numOverlap = entrezOverlap, numPresent = entrezPresent)
x = x[colorDat$col,]
x$numGenes = colorDat$numGenes
x$totGene = nrow(rse_gene)

## pink enrich
mat = matrix(c(5,198, 128, 13894 -5 - 198-128), nr=2)
chisq.test(mat)
getOR(mat)

magmaSig = magma[which(magma$P_JOINT
## line up
mm = match(rowData(rse_gene)$Symbol, magma$SYMBOL)
magGenes = split(magma$SYMBOL[mm], net$colorsLab)
