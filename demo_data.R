##
library(jaffelab)
library(SummarizedExperiment)
library(edgeR)

# num of reads, unfiltered
pd = read.csv("preprocessed_data/read_and_alignment_metrics_zandiHyde_Bipolar_LIBD.csv",
	as.is=TRUE,row.names=1)
quantile(pd$numReads/1e6)
	
# postQC 
load("data/zandiHypde_bipolar_rseGene_n511.rda")
getRPKM = recount::getRPKM
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx == "Other"] = "Bipolar"
pd$PrimaryDx = factor(as.character(pd$PrimaryDx), 
	levels = c("Control","Bipolar"))
# get RPKM
geneRpkm = getRPKM(rse_gene, length="Length")

# do PCA
pca = prcomp(t(log2(geneRpkm+1)))
pcaVars = getPcaVars(pca)

pdf("qc_checks/pcs_versus_stuff.pdf")
par(mar=c(11,6,2,2), cex.axis=2,cex.lab=2)
palette(brewer.pal(8,"Set1"))
## assignment PC1
plot(pca$x[,1] ~ pd$totalAssignedGene,
	xlab = "Gene Assignment Rate",
	pch = ifelse(pd$BrainRegion == "sACC", 21, 24),
	bg = factor(pd$PrimaryDx),
	ylab = paste0("PC1: ", pcaVars[1], "% Var Expl"))
legend("topright", levels(pd$PrimaryDx), col = 1:2, pch=15,cex=2)
legend("bottomleft", c("sACC","Amygdala"), col = "black",
	pt.bg="black", pch=c(21,24),cex=1.7)

## assignment PC1
plot(pca$x[,2] ~ pd$totalAssignedGene,
	xlab = "Gene Assignment Rate",
	pch = ifelse(pd$BrainRegion == "sACC", 21, 24),
	bg = factor(pd$PrimaryDx),
	ylab = paste0("PC2: ", pcaVars[2], "% Var Expl"))
legend("topright", levels(pd$PrimaryDx), col = 1:2, pch=15,cex=2)
legend("bottomleft", c("sACC","Amygdala"), col = "black",
	pt.bg="black", pch=c(21,24),cex=1.7)
		
## ercc
boxplot(pd$ERCCsumLogErr ~ ss(pd$SAMPLE_ID,"_",2),
	xlab = "",	ylab = "ERCC Bias",las=3,outline=FALSE)
points(pd$ERCCsumLogErr ~ jitter(as.numeric(factor(ss(pd$SAMPLE_ID,"_",2))),
	amount=0.15), pch = ifelse(pd$BrainRegion == "sACC", 21, 24),
	bg = factor(pd$PrimaryDx))
dev.off()

####
## control analysis for region
rse_gene_control = rse_gene[,rse_gene$PrimaryDx == "Control"]
rse_gene_control = rse_gene_control[rowMeans(getRPKM(rse_gene_control,"Length")) > 0.2,]

dge = DGEList(counts = assays(rse_gene_control)$counts, 
	genes = rowData(rse_gene_control))
dge = calcNormFactors(dge)

## mean-variance
mod0 = model.matrix(~BrainRegion, data=colData(rse_gene_control))
vGene0 = voom(dge,mod0,plot=TRUE)
fitGene0 = lmFit(vGene0)
ebGene0 = ebayes(fitGene0)

# w/ assignment adjusting
mod = model.matrix(~BrainRegion + totalAssignedGene,
	data=colData(rse_gene_control))
vGene = voom(dge,mod,plot=TRUE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)

plot(fitGene0$coef[,2], fitGene$coef[,2])
plot(ebGene0$t[,2], ebGene$t[,2])
plot(-log10(ebGene0$p[,2]), -log10(ebGene$p[,2]))

sum(p.adjust(ebGene$p[,2], "bonf") < 0.05)
sum(p.adjust(ebGene0$p[,2], "bonf") < 0.05)

pdf("qc_checks/volcano_region_controls_adjusted.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(fitGene$coef[,2], -log10(ebGene$p[,2]),
	pch = 21, bg=ifelse(p.adjust(ebGene$p[,2],
		"bonf") < 0.05, "red","grey"),
	xlab = "log2FC (sACC vs Amygdala)",
	ylab = "-log10(p-value)")
dev.off()
## do analysis
ebGene = ebayes(fitGene)

##################
# look at eQTLs  #
##################

load("eqtl_tables/mergedEqtl_output_gene_sacc.rda")
geneEqtl_sACC = geneEqtl
load("eqtl_tables/mergedEqtl_output_gene_amygdala.rda")
geneEqtl_Amygdala = geneEqtl
rm(geneEqtl)

## together
sigList = list(sACC = geneEqtl_sACC[geneEqtl_sACC$FDR < 0.01,],
	Amygdala = geneEqtl_Amygdala[geneEqtl_Amygdala$FDR < 0.01,])

sapply(sigList,nrow)

sapply(sigList,function(x) length(unique(x$snps)))
sapply(sigList,function(x) length(unique(x$gene)))
sapply(sigList,function(x) length(unique(x$Symbol)))

sapply(sigList,function(x) quantile(abs(x$beta)))

## pgc examples
geneEqtl_Amygdala[which(geneEqtl_Amygdala$snps == "rs7296288"),]
geneEqtl_Amygdala[which(geneEqtl_Amygdala$snps == "rs7296288"),]
geneEqtl_sACC[which(geneEqtl_sACC$snps == "rs12576775"),]

## load genotypes
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
s = snp["rs7296288",pd$BrNum]
mds = mds[pd$BrNum,1:5]

s[s==0] = "AA"
s[s==1] = "AC"
s[s==2] = "CC"
geneExprs = log2(geneRpkm+1)


pdf("qc_checks/example_eqtl.pdf",w=10,h=6)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
boxplot(geneExprs["ENSG00000167553.14",] ~ s*pd$BrainRegion,
	names = rep(sort(unique(s)),2),
	xlab = "rs7296288", ylab = "TUBA1C log2 Exprs")
text(x=c(2,5),y=4.6,c("Amygdala","sACC"),cex=2)
abline(v=3.5,lty=2)
boxplot(geneExprs["ENSG00000167550.10",] ~ s*pd$BrainRegion,
	names = rep(sort(unique(s)),2),
	xlab = "rs7296288", ylab = "RHEBL1 log2 Exprs")
abline(v=3.5,lty=2)
text(x=c(2,5),y=2.6,c("Amygdala","sACC"),cex=2)
dev.off()

## quick conditional

## overall
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	pd$BrainRegion+ pca$x[,1:10] + as.matrix(mds)))$coef[2,]
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	pd$BrainRegion +geneExprs["ENSG00000167550.10",] + pca$x[,1:10] + as.matrix(mds)))$coef[2,]	

summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	pd$BrainRegion+ pca$x[,1:10] + as.matrix(mds)))$coef[2,]
summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	pd$BrainRegion +geneExprs["ENSG00000167553.14",] + pca$x[,1:10] + as.matrix(mds)))$coef[2,]	
	
	
## sACC
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	pca$x[,1:10] + as.matrix(mds), subset=pd$BrainRegion=="sACC"))$coef[2,]
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	geneExprs["ENSG00000167550.10",] + pca$x[,1:10] + as.matrix(mds), 
		subset=pd$BrainRegion=="sACC"))$coef[2,]	

summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	pca$x[,1:10] + as.matrix(mds), subset=pd$BrainRegion=="sACC"))$coef[2,]
summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	geneExprs["ENSG00000167553.14",] + pca$x[,1:10] + as.matrix(mds), 
	subset=pd$BrainRegion=="sACC"))$coef[2,]	
	
		
## amygdala
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	pca$x[,1:10] + as.matrix(mds), subset=pd$BrainRegion=="Amygdala"))$coef[2,]
summary(lm(geneExprs["ENSG00000167553.14",] ~ as.numeric(factor(s)) + 
	geneExprs["ENSG00000167550.10",] + pca$x[,1:10] + as.matrix(mds), 
		subset=pd$BrainRegion=="Amygdala"))$coef[2,]	

summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	pca$x[,1:10] + as.matrix(mds), subset=pd$BrainRegion=="Amygdala"))$coef[2,]
summary(lm(geneExprs["ENSG00000167550.10",] ~ as.numeric(factor(s)) + 
	geneExprs["ENSG00000167553.14",] + pca$x[,1:10] + as.matrix(mds), 
	subset=pd$BrainRegion=="Amygdala"))$coef[2,]	
	
