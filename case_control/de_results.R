####
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(recount)
library(jaffelab)
library(IRanges)


## feature maps
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
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


################
## load

load("bipolarControl_deStats_byRegion_qSVAjoint.rda", verbose=TRUE)

## confirm order
feat = c(rownames(rse_gene), rownames(rse_exon),rownames(rse_jxn),rownames(rse_tx))
identical(rownames(statOut), feat)

## add featmap
statOut$Type = c(rep("1Gene",length(rse_gene)),rep("3Exon",length(rse_exon)),
				rep("4Jxn",length(rse_jxn)),rep("2Tx",length(rse_tx)))
statOut$gene_id = c(rowData(rse_gene)$ensemblID, rowData(rse_exon)$ensemblID,
					rowData(rse_jxn)$newGeneID, rowData(rse_tx)$gene_id)

					
## sig results
sigStat_amyg = statOut[which(statOut$"adj.P.Val_Amyg" < 0.05),]
sigStat_sacc = statOut[which(statOut$"adj.P.Val_sACC" < 0.05),]


################
## metrics

## total features
nrow(sigStat_amyg)  ## 3789
nrow(sigStat_sacc)   ## 1842

## per feature
table(sigStat_amyg$Type)
# 1Gene   2Tx 3Exon  4Jxn
# 104 	 3472   206     7

table(sigStat_sacc$Type)
# 1Gene   2Tx 3Exon  4Jxn
# 666   106   959   111


## unique ensemblIDs
length(unique(sigStat_amyg$gene_id))   ## 2903
length(unique(sigStat_sacc$gene_id))   ## 939

tapply(sigStat_amyg$gene_id, sigStat_amyg$Type, function(x) length(unique(x)))
# 1Gene   2Tx 3Exon  4Jxn
# 104  	 2756    70     6
tapply(sigStat_sacc$gene_id, sigStat_sacc$Type, function(x) length(unique(x)))
# 1Gene   2Tx 3Exon  4Jxn
# 666      90   388    56



gsigStat_amyg = sigStat_amyg[which(sigStat_amyg$Type=="Gene"),]
gsigStat_sacc = sigStat_sacc[which(sigStat_sacc$Type=="Gene"),]


library(VennDiagram)

venn.diagram(list(Amygdala = unique(gsigStat_amyg$gene_id), sACC = unique(gsigStat_sacc$gene_id)), 
	fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_qSVAjoint_gene.png")




# ## top gene
# gsigStat_amyg = sigStat_amyg[which(sigStat_amyg$Type=="1Gene"),]
# gsigStat_sacc = sigStat_sacc[which(sigStat_sacc$Type=="1Gene"),]

# pd = colData(rse_gene)
# pd_amyg = pd[which(pd$BrainRegion=="Amygdala"),]
# pd_sacc = pd[which(pd$BrainRegion=="sACC"),]

# geneRpkm = getRPKM(rse_gene)
# geneRpkm_amyg = geneRpkm[,which(pd$BrainRegion=="Amygdala")]
# geneRpkm_sacc = geneRpkm[,which(pd$BrainRegion=="sACC")]

# gsigStat_amyg = gsigStat_amyg[order(gsigStat_amyg$P.Value_Amyg, decreasing=FALSE),]
# gsigStat_sacc = gsigStat_sacc[order(gsigStat_sacc$P.Value_sACC, decreasing=FALSE),]


# topamyg = rownames(gsigStat_amyg)[1]

# pdf(paste0("topGeneDE.pdf"),h=6,w=8)
# par(mar=c(5,4,3,2),cex.axis=1.3,cex.lab=1.5,cex.main=1.5)
  # boxplot(geneRpkm_amyg[topamyg,] ~ pd_amyg$PrimaryDx, xlab="Dx", 
		# main=topamyg, outline=FALSE)
  # points(geneRpkm_amyg[topamyg,] ~ jitter(as.numeric(factor(pd_amyg$PrimaryDx)), 
		  # amount=0.15), cex=1.5,
		  # pch=21, bg = "blue")
# dev.off()

















