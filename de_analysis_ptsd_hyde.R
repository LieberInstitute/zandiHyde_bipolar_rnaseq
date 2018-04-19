######

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(sva)
library(limma)
library(edgeR)

## load DG
load("data/zandiHypde_bipolar_rseTx_n511.rda")
load("data/zandiHypde_bipolar_rseJxn_n511.rda")
load("data/zandiHypde_bipolar_rseExon_n511.rda")
load("data/zandiHypde_bipolar_rseGene_n511.rda")

## load SNP data for MDS
load("genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")

#### function for RPKM, in latest recount...
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colSums(assays(rse)$counts)), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}
getRPM = function(rse, target = 10e6) {
	require(SummarizedExperiment)
	bg = matrix(rep(colSums(assays(rse)$counts/target)), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

## load matching info for PTSD samples
ptsd_matches = read.csv("ptsd_matches.csv", row.names=1)
ptsd_matches = ptsd_matches[!is.na(ptsd_matches$MatchIndex),]

## subset rses to ptsd dataset
ptsdInd = which(colData(rse_gene)$RNum %in% ptsd_matches$RNum)
identical(colData(rse_gene)$SAMPLE_ID, colData(rse_jxn)$SAMPLE_ID)  ## confirm consistent order
# cov_rse = cov_rse[,ptsdInd]
rse_gene = rse_gene[,ptsdInd]
rse_exon = rse_exon[,ptsdInd]
rse_jxn = rse_jxn[,ptsdInd]
rse_tx = rse_tx[,ptsdInd]

## add PTSD to Dx column
colData(rse_gene)$PrimaryDx = colData(rse_exon)$PrimaryDx = 
		colData(rse_jxn)$PrimaryDx = colData(rse_tx)$PrimaryDx = 
		ptsd_matches[rownames(colData(rse_gene)),"PrimaryDx"]
## reorder		
colData(rse_gene)$PrimaryDx = colData(rse_exon)$PrimaryDx = 
		colData(rse_jxn)$PrimaryDx = colData(rse_tx)$PrimaryDx = 
		factor(colData(rse_gene)$PrimaryDx, levels = c("Bipolar_PTSD", "Bipolar", "Control"))

		
## separate region analyses		
amyInd = which(colData(rse_gene)$BrainRegion == "Amygdala")
saccInd = which(colData(rse_gene)$BrainRegion == "sACC")
rse_gene_amy = rse_gene[,amyInd]
rse_gene_sacc = rse_gene[,saccInd]

rse_exon_amy = rse_exon[,amyInd]
rse_exon_sacc = rse_exon[,saccInd]

rse_jxn_amy = rse_jxn[,amyInd]
rse_jxn_sacc = rse_jxn[,saccInd]

mds1 = mds[colData(rse_exon_amy)$BrNum,]
mds2 = mds[colData(rse_gene_sacc)$BrNum,]



######################
#### Amygdala ######
######################

######################
# statistical model ##
mod0 = model.matrix(~PrimaryDx, data = colData(rse_gene_amy))
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds1[,1:5]),
				data = colData(rse_gene_amy))

#################
## modeling #####

### gene 
dge = DGEList(counts = assays(rse_gene_amy)$counts, genes = rowData(rse_gene_amy))
dge = calcNormFactors(dge)

## do analysis
vGene0 = voom(dge,mod0,plot=TRUE)
fitGene0 = lmFit(vGene0)
ebGene0 = ebayes(fitGene0)

#########################
## adjust for stuff
vGene_Amy = voom(dge,mod,plot=TRUE)
fitGene_Amy = lmFit(vGene_Amy)
ebGene_Amy = ebayes(fitGene_Amy)

# ## compare t-stats
# plot(ebGene_Amy$t[,2], ebGene0$t[,2])   # PTSD vs Bipolar
# abline(0,1, col="gray") 
# plot(ebGene_Amy$t[,3], ebGene0$t[,3])   # PTSD vs Control
# abline(0,1, col="gray") 

############ exons ############

######################
# statistical model ##
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds1[,1:5]),
				data = colData(rse_exon_amy))
#################
## modeling #####

dge = DGEList(counts = assays(rse_exon_amy)$counts, genes = rowData(rse_exon_amy))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vExon_Amy = voom(dge,mod,plot=TRUE)
fitExon_Amy = lmFit(vExon_Amy)
ebExon_Amy = ebayes(fitExon_Amy)


############ junctions ############

######################
# statistical model ##
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds1[,1:5]),
				data = colData(rse_jxn_amy))
#################
## modeling #####

dge = DGEList(counts = assays(rse_jxn_amy)$counts, genes = rowData(rse_jxn_amy))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vJxn_Amy = voom(dge,mod,plot=TRUE)
fitJxn_Amy = lmFit(vJxn_Amy)
ebJxn_Amy = ebayes(fitJxn_Amy)





######################
#### sACC ######
######################

######################
# statistical model ##
mod0 = model.matrix(~PrimaryDx, data = colData(rse_gene_sacc))
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds2[,1:5]),
				data = colData(rse_gene_sacc))

#################
## modeling #####

### gene 
dge = DGEList(counts = assays(rse_gene_sacc)$counts, genes = rowData(rse_gene_sacc))
dge = calcNormFactors(dge)

## do analysis
vGene0 = voom(dge,mod0,plot=TRUE)
fitGene0 = lmFit(vGene0)
ebGene0 = ebayes(fitGene0)

#########################
## adjust for stuff
vGene_Sacc = voom(dge,mod,plot=TRUE)
fitGene_Sacc = lmFit(vGene_Sacc)
ebGene_Sacc = ebayes(fitGene_Sacc)

# ## compare t-stats
# plot(ebGene_Sacc$t[,2], ebGene0$t[,2])   # PTSD vs Bipolar
# abline(0,1, col="gray") 
# plot(ebGene_Sacc$t[,3], ebGene0$t[,3])   # PTSD vs Control
# abline(0,1, col="gray") 


############ exons ############

######################
# statistical model ##
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds2[,1:5]),
				data = colData(rse_exon_sacc))
#################
## modeling #####

dge = DGEList(counts = assays(rse_exon_sacc)$counts, genes = rowData(rse_exon_sacc))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vExon_Sacc = voom(dge,mod,plot=TRUE)
fitExon_Sacc = lmFit(vExon_Sacc)
ebExon_Sacc = ebayes(fitExon_Sacc)


############ junctions ############

######################
# statistical model ##
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds2[,1:5]),
				data = colData(rse_jxn_sacc))
#################
## modeling #####

dge = DGEList(counts = assays(rse_jxn_sacc)$counts, genes = rowData(rse_jxn_sacc))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vJxn_Sacc = voom(dge,mod,plot=TRUE)
fitJxn_Sacc = lmFit(vJxn_Sacc)
ebJxn_Sacc = ebayes(fitJxn_Sacc)





identical(rownames(ebGene_Sacc$t),rownames(ebGene_Amy$t) )
identical(rownames(ebExon_Sacc$t),rownames(ebExon_Amy$t) )



######################
#### DG ######
######################

## load DG
load("/dcl01/ajaffe/data/lab/dg_hippo/count_data/astellas_dg_hg38_rseTx_n263.rda")
load("/dcl01/ajaffe/data/lab/dg_hippo/count_data/astellas_dg_hg38_rseJxn_n263.rda")
load("/dcl01/ajaffe/data/lab/dg_hippo/count_data/astellas_dg_hg38_rseExon_n263.rda")
load("/dcl01/ajaffe/data/lab/dg_hippo/count_data/astellas_dg_hg38_rseGene_n263.rda")

identical(rowData(rse_gene)$gencodeID, rowData(rse_gene_amy)$gencodeID)  ## same gene order as other regions?
identical(rowData(rse_exon)$gencodeID, rowData(rse_exon_amy)$gencodeID)

## load SNP data for MDS
load("/dcl01/ajaffe/data/lab/dg_hippo/genotype_data/astellas_dg_genotype_data_n263.rda")

## load matching info for PTSD samples
ptsd_matches = read.csv("/dcl01/ajaffe/data/lab/dg_hippo/ptsd_matches.csv", row.names=1)
ptsd_matches = ptsd_matches[!is.na(ptsd_matches$MatchIndex),]

## subset rses to ptsd dataset
ptsdInd = which(colData(rse_gene)$RNum %in% ptsd_matches$RNum)
rse_gene = rse_gene[,ptsdInd]
rse_exon = rse_exon[,ptsdInd]
rse_jxn = rse_jxn[,ptsdInd]
rse_tx = rse_tx[,ptsdInd]

## add PTSD to Dx column	
colData(rse_gene)$Dx = colData(rse_exon)$Dx = 
		colData(rse_jxn)$Dx = colData(rse_tx)$Dx = 
		ptsd_matches[rownames(colData(rse_gene)),"Dx"]
## reorder
colData(rse_gene)$Dx = colData(rse_exon)$Dx = 
		colData(rse_jxn)$Dx = colData(rse_tx)$Dx = 
		factor(colData(rse_gene)$Dx, levels = c("Bipolar_PTSD", "Bipolar", "Control"))
		
mds = mds[colData(rse_gene)$BrNum,]

######################
# statistical model ##
######################
mod0 = model.matrix(~Dx, data = colData(rse_gene))
mod = model.matrix(~Dx + totalAssignedGene + Sex + Age + 
				mitoRate + as.matrix(mds[,1:5]),
				data = colData(rse_gene))

#################
## modeling #####

### gene 
dge = DGEList(counts = assays(rse_gene)$counts, genes = rowData(rse_gene))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vGene_DG = voom(dge,mod,plot=TRUE)
fitGene_DG = lmFit(vGene_DG)
ebGene_DG = ebayes(fitGene_DG)



############ exons ############
mod = model.matrix(~Dx + totalAssignedGene + Sex + Age + 
				mitoRate + as.matrix(mds[,1:5]),
				data = colData(rse_exon))
#################
## modeling #####

dge = DGEList(counts = assays(rse_exon)$counts, genes = rowData(rse_exon))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vExon_DG = voom(dge,mod,plot=TRUE)
fitExon_DG = lmFit(vExon_DG)
ebExon_DG = ebayes(fitExon_DG)


############ junctions ############
mod = model.matrix(~Dx + totalAssignedGene + Sex + Age + 
				mitoRate + as.matrix(mds[,1:5]),
				data = colData(rse_jxn))
#################
## modeling #####

dge = DGEList(counts = assays(rse_jxn)$counts, genes = rowData(rse_jxn))
dge = calcNormFactors(dge)

#########################
## adjust for stuff
vJxn_DG = voom(dge,mod,plot=TRUE)
fitJxn_DG = lmFit(vJxn_DG)
ebJxn_DG = ebayes(fitJxn_DG)



## one last order check
identical(rownames(ebGene_DG$t),rownames(ebGene_Sacc$t))
identical(rownames(ebExon_DG$t),rownames(ebExon_Sacc$t))
identical(rownames(ebJxn_DG$t),rownames(ebJxn_Sacc$t)) ## FALSE - different junctions



###########################
# combine region metrics ##
###########################

##### metrics ######
geneStats = as.data.frame(rowRanges(rse_gene))
names(geneStats)[1:5] = c("Chr","Start","End","Width","Strand")
geneStats$Length = geneStats$Class = geneStats$meanExprs = geneStats$gencodeTx = NULL

geneStats$Amy_meanRpkm = rowMeans(getRPKM(rse_gene_amy))
geneStats$Sacc_meanRpkm = rowMeans(getRPKM(rse_gene_sacc))
geneStats$DG_meanRpkm = rowMeans(getRPKM(rse_gene))


##### Amygdala
## fold change
geneStats$Amy_log2FC_PTSD_BP = fitGene_Amy$coef[,2]
geneStats$Amy_log2FC_PTSD_CNT = fitGene_Amy$coef[,3]
## tstat 
geneStats$Amy_tstat_PTSD_BP = ebGene_Amy$t[,2]
geneStats$Amy_tstat_PTSD_CNT = ebGene_Amy$t[,3]
## pvalue 
geneStats$Amy_pvalue_PTSD_BP = ebGene_Amy$p[,2]
geneStats$Amy_pvalue_PTSD_CNT = ebGene_Amy$p[,3]
## qvalue 
geneStats$Amy_qvalue_PTSD_BP = NA
geneStats$Amy_qvalue_PTSD_BP[geneStats$Amy_meanRpkm>0] = p.adjust(geneStats$Amy_pvalue_PTSD_BP[geneStats$Amy_meanRpkm>0] , "fdr")
geneStats$Amy_qvalue_PTSD_CNT = NA
geneStats$Amy_qvalue_PTSD_CNT[geneStats$Amy_meanRpkm>0] = p.adjust(geneStats$Amy_pvalue_PTSD_CNT[geneStats$Amy_meanRpkm>0] , "fdr")


##### sACC
## fold change
geneStats$Sacc_log2FC_PTSD_BP = fitGene_Sacc$coef[,2]
geneStats$Sacc_log2FC_PTSD_CNT = fitGene_Sacc$coef[,3]
## tstat 
geneStats$Sacc_tstat_PTSD_BP = ebGene_Sacc$t[,2]
geneStats$Sacc_tstat_PTSD_CNT = ebGene_Sacc$t[,3]
## pvalue 
geneStats$Sacc_pvalue_PTSD_BP = ebGene_Sacc$p[,2]
geneStats$Sacc_pvalue_PTSD_CNT = ebGene_Sacc$p[,3]
## qvalue 
geneStats$Sacc_qvalue_PTSD_BP = NA
geneStats$Sacc_qvalue_PTSD_BP[geneStats$Sacc_meanRpkm>0] = p.adjust(geneStats$Sacc_pvalue_PTSD_BP[geneStats$Sacc_meanRpkm>0] , "fdr")
geneStats$Sacc_qvalue_PTSD_CNT = NA
geneStats$Sacc_qvalue_PTSD_CNT[geneStats$Sacc_meanRpkm>0] = p.adjust(geneStats$Sacc_pvalue_PTSD_CNT[geneStats$Sacc_meanRpkm>0] , "fdr")


##### DG
## fold change
geneStats$DG_log2FC_PTSD_BP = fitGene_DG$coef[,2]
geneStats$DG_log2FC_PTSD_CNT = fitGene_DG$coef[,3]
## tstat 
geneStats$DG_tstat_PTSD_BP = ebGene_DG$t[,2]
geneStats$DG_tstat_PTSD_CNT = ebGene_DG$t[,3]
## pvalue 
geneStats$DG_pvalue_PTSD_BP = ebGene_DG$p[,2]
geneStats$DG_pvalue_PTSD_CNT = ebGene_DG$p[,3]
## qvalue 
geneStats$DG_qvalue_PTSD_BP = NA
geneStats$DG_qvalue_PTSD_BP[geneStats$DG_meanRpkm>0] = p.adjust(geneStats$DG_pvalue_PTSD_BP[geneStats$DG_meanRpkm>0] , "fdr")
geneStats$DG_qvalue_PTSD_CNT = NA
geneStats$DG_qvalue_PTSD_CNT[geneStats$DG_meanRpkm>0] = p.adjust(geneStats$DG_pvalue_PTSD_CNT[geneStats$DG_meanRpkm>0] , "fdr")


write.csv(geneStats, "ptsd_geneStats2.csv")





###########################
######### boxplots ########
###########################
library(RColorBrewer)

### genes to plot		
genes = read.csv("ptsd_geneStats_filtered.csv", row.names=1 )

## residualize expression - Amygdala		
amyRpkm = getRPKM(rse_gene_amy)
yExprs_amy = log2(amyRpkm+1)
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds1[,1:5]),
				data = colData(rse_gene_amy))				
yExprsAdj_amy = cleaningY(yExprs_amy, mod, P=3)

## residualize expression - sACC	
saccRpkm = getRPKM(rse_gene_sacc)
yExprs_sacc = log2(saccRpkm+1)
mod = model.matrix(~PrimaryDx + totalAssignedGene + Sex + AgeDeath + 
				mitoRate + as.matrix(mds2[,1:5]),
				data = colData(rse_gene_sacc))			
yExprsAdj_sacc = cleaningY(yExprs_sacc, mod, P=3)

## residualize expression - DG	
dgRpkm = getRPKM(rse_gene)
yExprs_dg = log2(dgRpkm+1)
mod = model.matrix(~Dx + totalAssignedGene + Sex + Age + 
				mitoRate + as.matrix(mds[,1:5]),
				data = colData(rse_gene))				
yExprsAdj_dg = cleaningY(yExprs_dg, mod, P=3)


yExprs = cbind(yExprsAdj_amy,yExprsAdj_sacc,yExprsAdj_dg)
yExprs = yExprs[rownames(genes),]   ## only plots genes from csv list


### combine colData from regions
pd = data.frame(samp = colnames(yExprs),
			region = c(as.character(colData(rse_gene_amy)$BrainRegion), 
						as.character(colData(rse_gene_sacc)$BrainRegion), 
						as.character(colData(rse_gene)$Region)),
			dx = c(as.character(colData(rse_gene_amy)$PrimaryDx), 
						as.character(colData(rse_gene_sacc)$PrimaryDx), 
						as.character(colData(rse_gene)$Dx))
				)
pd$region = as.character(pd$region)
pd$region[pd$region=="DentateGyrus"] = "DG"

pd$dx = as.character(pd$dx)
pd$dx[pd$dx=="Bipolar_PTSD"] = "PTSD"

pd$group = paste0(pd$dx, "_", pd$region)
pd$group = factor(pd$group, 
		levels=c("Control_Amygdala","PTSD_Amygdala","Bipolar_Amygdala",
		"Control_sACC","PTSD_sACC","Bipolar_sACC",
		"Control_DG","PTSD_DG","Bipolar_DG" ) )


	

### make plots
pdf("ptsd_geneStats_filtered_adj.pdf", h=6, w=10)
par(mar=c(7,6,4,2),cex.axis=1,cex.lab=1.2,cex.main=1.5)
palette(brewer.pal(6,"Spectral"))
				
for (i in 1:nrow(yExprs)) {

boxplot(yExprs[i,] ~ pd$group,
		xlab="", main=paste0(genes$gencodeID[i],"\n",genes$Symbol[i]),
		ylab="Residualized Expression", las=2, outline=FALSE,
		xaxt="n")
points(yExprs[i,] ~ jitter(as.numeric(factor(pd$group))),
		   pch=21, 
		   bg=as.numeric(as.factor(pd$region))*2,cex=1.8)
abline(v=c(3.5,6.5),lty=2, col="gray")

axis(1, at=1:9, labels=rep(c("Control","PTSD","Bipolar"),3) )
axis(1, at=c(2,5,8), labels=c("Amygdala","sACC","DG"), line=2.5, tick=FALSE, cex.axis=1.2 )

legend("topleft", paste0("p=",signif(genes$Amy_pvalue_PTSD_CNT[i],3)) )
legend("top", paste0("p=",signif(genes$Sacc_pvalue_PTSD_CNT[i],3)) )
legend("topright", paste0("p=",signif(genes$DG_pvalue_PTSD_CNT[i],3)) )
}

dev.off()










#################
##### exons #####
#################

##### metrics ######
exonStats = as.data.frame(rowRanges(rse_exon))
names(exonStats)[1:5] = c("Chr","Start","End","Width","Strand")
exonStats$Length = exonStats$Class = exonStats$meanExprs = exonStats$gencodeTx = NULL

exonStats$Amy_meanRpkm = rowMeans(getRPKM(rse_exon_amy))
exonStats$Sacc_meanRpkm = rowMeans(getRPKM(rse_exon_sacc))
exonStats$DG_meanRpkm = rowMeans(getRPKM(rse_exon))


##### Amygdala
## fold change
exonStats$Amy_log2FC_PTSD_BP = fitExon_Amy$coef[,2]
exonStats$Amy_log2FC_PTSD_CNT = fitExon_Amy$coef[,3]
## tstat 
exonStats$Amy_tstat_PTSD_BP = ebExon_Amy$t[,2]
exonStats$Amy_tstat_PTSD_CNT = ebExon_Amy$t[,3]
## pvalue 
exonStats$Amy_pvalue_PTSD_BP = ebExon_Amy$p[,2]
exonStats$Amy_pvalue_PTSD_CNT = ebExon_Amy$p[,3]
## qvalue 
exonStats$Amy_qvalue_PTSD_BP = NA
exonStats$Amy_qvalue_PTSD_BP[exonStats$Amy_meanRpkm>0] = p.adjust(exonStats$Amy_pvalue_PTSD_BP[exonStats$Amy_meanRpkm>0] , "fdr")
exonStats$Amy_qvalue_PTSD_CNT = NA
exonStats$Amy_qvalue_PTSD_CNT[exonStats$Amy_meanRpkm>0] = p.adjust(exonStats$Amy_pvalue_PTSD_CNT[exonStats$Amy_meanRpkm>0] , "fdr")


##### sACC
## fold change
exonStats$Sacc_log2FC_PTSD_BP = fitExon_Sacc$coef[,2]
exonStats$Sacc_log2FC_PTSD_CNT = fitExon_Sacc$coef[,3]
## tstat 
exonStats$Sacc_tstat_PTSD_BP = ebExon_Sacc$t[,2]
exonStats$Sacc_tstat_PTSD_CNT = ebExon_Sacc$t[,3]
## pvalue 
exonStats$Sacc_pvalue_PTSD_BP = ebExon_Sacc$p[,2]
exonStats$Sacc_pvalue_PTSD_CNT = ebExon_Sacc$p[,3]
## qvalue 
exonStats$Sacc_qvalue_PTSD_BP = NA
exonStats$Sacc_qvalue_PTSD_BP[exonStats$Sacc_meanRpkm>0] = p.adjust(exonStats$Sacc_pvalue_PTSD_BP[exonStats$Sacc_meanRpkm>0] , "fdr")
exonStats$Sacc_qvalue_PTSD_CNT = NA
exonStats$Sacc_qvalue_PTSD_CNT[exonStats$Sacc_meanRpkm>0] = p.adjust(exonStats$Sacc_pvalue_PTSD_CNT[exonStats$Sacc_meanRpkm>0] , "fdr")


##### DG
## fold change
exonStats$DG_log2FC_PTSD_BP = fitExon_DG$coef[,2]
exonStats$DG_log2FC_PTSD_CNT = fitExon_DG$coef[,3]
## tstat 
exonStats$DG_tstat_PTSD_BP = ebExon_DG$t[,2]
exonStats$DG_tstat_PTSD_CNT = ebExon_DG$t[,3]
## pvalue 
exonStats$DG_pvalue_PTSD_BP = ebExon_DG$p[,2]
exonStats$DG_pvalue_PTSD_CNT = ebExon_DG$p[,3]
## qvalue 
exonStats$DG_qvalue_PTSD_BP = NA
exonStats$DG_qvalue_PTSD_BP[exonStats$DG_meanRpkm>0] = p.adjust(exonStats$DG_pvalue_PTSD_BP[exonStats$DG_meanRpkm>0] , "fdr")
exonStats$DG_qvalue_PTSD_CNT = NA
exonStats$DG_qvalue_PTSD_CNT[exonStats$DG_meanRpkm>0] = p.adjust(exonStats$DG_pvalue_PTSD_CNT[exonStats$DG_meanRpkm>0] , "fdr")


write.csv(geneStats, "ptsd_geneStats2.csv")



################
##### jxns #####
#################

##### metrics ######
jxnStatsDG = as.data.frame(rowRanges(rse_jxn))
names(jxnStatsDG)[1:5] = c("Chr","Start","End","Width","Strand")
jxnStatsDG$meanExprs = jxnStatsDG$gencodeTx = NULL
jxnStatsDG$DG_meanRpkm = rowMeans(getRPM(rse_jxn))

jxnStatsSA = as.data.frame(rowRanges(rse_jxn_amy))
names(jxnStatsSA)[1:5] = c("Chr","Start","End","Width","Strand")
jxnStatsSA$meanExprs = jxnStatsSA$gencodeTx = NULL
jxnStatsSA$Amy_meanRpkm = rowMeans(getRPM(rse_jxn_amy))
jxnStatsSA$Sacc_meanRpkm = rowMeans(getRPM(rse_jxn_sacc))


##### Amygdala
## fold change
jxnStatsSA$Amy_log2FC_PTSD_BP = fitJxn_Amy$coef[,2]
jxnStatsSA$Amy_log2FC_PTSD_CNT = fitJxn_Amy$coef[,3]
## tstat 
jxnStatsSA$Amy_tstat_PTSD_BP = ebJxn_Amy$t[,2]
jxnStatsSA$Amy_tstat_PTSD_CNT = ebJxn_Amy$t[,3]
## pvalue 
jxnStatsSA$Amy_pvalue_PTSD_BP = ebJxn_Amy$p[,2]
jxnStatsSA$Amy_pvalue_PTSD_CNT = ebJxn_Amy$p[,3]
## qvalue 
jxnStatsSA$Amy_qvalue_PTSD_BP = NA
jxnStatsSA$Amy_qvalue_PTSD_BP[jxnStatsSA$Amy_meanRpkm>0] = p.adjust(jxnStatsSA$Amy_pvalue_PTSD_BP[jxnStatsSA$Amy_meanRpkm>0] , "fdr")
jxnStatsSA$Amy_qvalue_PTSD_CNT = NA
jxnStatsSA$Amy_qvalue_PTSD_CNT[jxnStatsSA$Amy_meanRpkm>0] = p.adjust(jxnStatsSA$Amy_pvalue_PTSD_CNT[jxnStatsSA$Amy_meanRpkm>0] , "fdr")


##### sACC
## fold change
jxnStatsSA$Sacc_log2FC_PTSD_BP = fitJxn_Sacc$coef[,2]
jxnStatsSA$Sacc_log2FC_PTSD_CNT = fitJxn_Sacc$coef[,3]
## tstat 
jxnStatsSA$Sacc_tstat_PTSD_BP = ebJxn_Sacc$t[,2]
jxnStatsSA$Sacc_tstat_PTSD_CNT = ebJxn_Sacc$t[,3]
## pvalue 
jxnStatsSA$Sacc_pvalue_PTSD_BP = ebJxn_Sacc$p[,2]
jxnStatsSA$Sacc_pvalue_PTSD_CNT = ebJxn_Sacc$p[,3]
## qvalue 
jxnStatsSA$Sacc_qvalue_PTSD_BP = NA
jxnStatsSA$Sacc_qvalue_PTSD_BP[jxnStatsSA$Sacc_meanRpkm>0] = p.adjust(jxnStatsSA$Sacc_pvalue_PTSD_BP[jxnStatsSA$Sacc_meanRpkm>0] , "fdr")
jxnStatsSA$Sacc_qvalue_PTSD_CNT = NA
jxnStatsSA$Sacc_qvalue_PTSD_CNT[jxnStatsSA$Sacc_meanRpkm>0] = p.adjust(jxnStatsSA$Sacc_pvalue_PTSD_CNT[jxnStatsSA$Sacc_meanRpkm>0] , "fdr")


##### DG
## fold change
jxnStatsDG$DG_log2FC_PTSD_BP = fitJxn_DG$coef[,2]
jxnStatsDG$DG_log2FC_PTSD_CNT = fitJxn_DG$coef[,3]
## tstat 
jxnStatsDG$DG_tstat_PTSD_BP = ebJxn_DG$t[,2]
jxnStatsDG$DG_tstat_PTSD_CNT = ebJxn_DG$t[,3]
## pvalue 
jxnStatsDG$DG_pvalue_PTSD_BP = ebJxn_DG$p[,2]
jxnStatsDG$DG_pvalue_PTSD_CNT = ebJxn_DG$p[,3]
## qvalue 
jxnStatsDG$DG_qvalue_PTSD_BP = NA
jxnStatsDG$DG_qvalue_PTSD_BP[jxnStatsDG$DG_meanRpkm>0] = p.adjust(jxnStatsDG$DG_pvalue_PTSD_BP[jxnStatsDG$DG_meanRpkm>0] , "fdr")
jxnStatsDG$DG_qvalue_PTSD_CNT = NA
jxnStatsDG$DG_qvalue_PTSD_CNT[jxnStatsDG$DG_meanRpkm>0] = p.adjust(jxnStatsDG$DG_pvalue_PTSD_CNT[jxnStatsDG$DG_meanRpkm>0] , "fdr")


geneB = geneStats[which(geneStats$Symbol=="BDNF"),]
exonB = exonStats[which(exonStats$Symbol=="BDNF"),]
geneexon = rbind(geneB, exonB)

jxnB1 = jxnStatsSA[which(jxnStatsSA$Symbol=="BDNF"),]
jxnB2 = jxnStatsDG[which(jxnStatsDG$Symbol=="BDNF"),]


write.csv(geneexon, "ptsd_bdnf_gene_exon.csv")
write.csv(jxnB1, "ptsd_bdnf_jxn1.csv")
write.csv(jxnB2, "ptsd_bdnf_jxn2.csv")





