
###############################################################################################################################################################################
##differentially expressed features: genes, exons, junctions, transcripts
#/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda

#I used gene symbols for the following analyses. Note, this excludes some features with ensemblIDs, etc, but no gene symbol.

##Note, CSEA may not be ideal for features other than genes because each gene symbol is only counted once
#if multiple transcripts map to the same gene, perhaps those additional hits should be factored in
###############################################################################################################################################################################

####The significance values span wide ranges for the different features as well as the different regions
###See bipolar_csea_code.R for sensitivity analyses evaluating significance thresholding and the effect on cell type results
##Results of sensitivity analyses are in /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses
#Thus, I established cutoffs for pSI input based on the number of gene symbols, not a significance threshold

library(jaffelab)
library(SummarizedExperiment)
library(pSI)
setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/")

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda", verbose=TRUE)

#Only keep class = InGen and symbol != blank so we don't have issues later with selecting exactly 225 genes
statOut <- statOut[statOut$Class == "InGen" & statOut$Symbol != "",]

#Subset the features
statOutGene <- statOut[statOut$Type == "Gene",]
statOutExon <- statOut[statOut$Type == "Exon",]
statOutJxn <- statOut[statOut$Type == "Junction",]
statOutTx <- statOut[statOut$Type == "Transcript",]
dim(statOutGene)
#[1] 19205    17
dim(statOutExon)
#[1] 382411     17
dim(statOutJxn)
#[1] 208163     17
dim(statOutTx) 
#[1] 73214    17


#pSI calculations work well with about 200-250 genes, based on experience and the original publication (doi: 10.1523/JNEUROSCI.4488-13.2014)
#though you do have to consider, when thinking about pSIs, that both input gene number and input gene purity (for a given cell type) will affect results
#Here, we'll take the top 225 unique genes based on ordered q values
#Note, some of these features will have a gencodeID, etc. but no gene symbol, and since CSEA takes gene symbols only, I took the top 225 unique gene symbol names


#Order objects based on q values
statOutGeneAmyg <- statOutGene[order(statOutGene$adj.P.Val_Amyg),]
statOutExonAmyg <- statOutExon[order(statOutExon$adj.P.Val_Amyg),]
statOutJxnAmyg <- statOutJxn[order(statOutJxn$adj.P.Val_Amyg),]
statOutTxAmyg <- statOutTx[order(statOutTx$adj.P.Val_Amyg),]

statOutGenesACC <- statOutGene[order(statOutGene$adj.P.Val_sACC),]
statOutExonsACC <- statOutExon[order(statOutExon$adj.P.Val_sACC),]
statOutJxnsACC <- statOutJxn[order(statOutJxn$adj.P.Val_sACC),]
statOutTxsACC <- statOutTx[order(statOutTx$adj.P.Val_sACC),]

#find number of rows needed for 225 unique genes 
length(unique(statOutGeneAmyg[c(1:225),"Symbol"]))
length(unique(statOutExonAmyg[c(1:577),"Symbol"]))
length(unique(statOutJxnAmyg[c(1:312),"Symbol"]))
length(unique(statOutTxAmyg[c(1:239),"Symbol"]))

length(unique(statOutGenesACC[c(1:225),"Symbol"]))
length(unique(statOutExonsACC[c(1:611),"Symbol"]))
length(unique(statOutJxnsACC[c(1:362),"Symbol"]))
length(unique(statOutTxsACC[c(1:237),"Symbol"]))


#Name objects
gene_amyg_225 <- statOutGeneAmyg[c(1:225), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
exon_amyg_225 <- statOutExonAmyg[c(1:577), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_amyg_225 <- statOutJxnAmyg[c(1:312), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
tx_amyg_225 <- statOutTxAmyg[c(1:239), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]

gene_sacc_225 <- statOutGenesACC[c(1:225), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
exon_sacc_225 <- statOutExonsACC[c(1:611), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_sacc_225 <- statOutJxnsACC[c(1:362), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
tx_sacc_225 <- statOutTxsACC[c(1:237), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]

#Save objects as rdas and csvs
save(gene_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
save(exon_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
save(jxn_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
save(tx_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")

save(gene_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
save(exon_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
save(jxn_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
save(tx_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")

write.csv(gene_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_225.csv")
write.csv(exon_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_Amyg_225.csv")
write.csv(jxn_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_Amyg_225.csv")
write.csv(tx_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_Amyg_225.csv")

write.csv(gene_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_225.csv")
write.csv(exon_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_sACC_225.csv")
write.csv(jxn_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_sACC_225.csv")
write.csv(tx_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_sACC_225.csv")



##### Xu and Dougherty et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### CSEA
### Online tool (http://genetics.wustl.edu/jdlab/csea-tool-2/)
##Use rdas generated above

#CSEA using the pSI R package, built from source (http://genetics.wustl.edu/jdlab/psi_package/)
library(gdata)
library(pSI)
library(pSI.data)
data(mouse)
data(human)
library(tibble)
library(dplyr)

#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda",verbose=TRUE)
gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(mouse$psi.out,gene_amyg_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(mouse$psi.out,exon_amyg_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,jxn_amyg_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,tx_amyg_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	

write.csv(fisher.iteration(mouse$psi.out,gene_amyg_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,exon_amyg_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,jxn_amyg_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,tx_amyg_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Dougherty.csv")



#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda",verbose=TRUE)
gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(mouse$psi.out,gene_sacc_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
			
fisher.iteration(mouse$psi.out,exon_sacc_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(mouse$psi.out,jxn_sacc_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(mouse$psi.out,tx_sacc_225_sym,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
write.csv(fisher.iteration(mouse$psi.out,gene_sacc_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,exon_sacc_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,jxn_sacc_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,tx_sacc_225_sym,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Dougherty.csv")		


##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.

#get data proper way
#process raw matrix
library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(SingleCellExperiment)
library(scater)
library(tibble)
library(dplyr)

setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/")

#################################################################################################################################
#get data for pSI calculation
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Sestan.adultHumanNuclei.Psychencode.Rdata",verbose=TRUE)
#the meta2 file in the above indicates which nuclei are left after QC filtering. ill just match those nuclei numbers to the raw counts matrix to subset. 
#note, Lake 2018's "raw" matrix file is already subsetted to those nuclei that have passed QC

umi.raw.subset <- umi.raw[,which(colnames(umi.raw) %in% rownames(meta2))]
#get rownames correctly formatted as well
rownames(umi.raw.subset) <- sub(".*\\|","",rownames(umi.raw.subset))

#drop any genes with no expression across all cells/nuclei
min(umi.raw.subset)
sum(is.na(umi.raw.subset))
totals <- rowSums(umi.raw.subset)
min(totals)
#appears no genes have zero counts across all nuclei (any genes that did were probably dropped in Li 2018's QC), but do below to check
umi.raw.subset_notzero <- umi.raw.subset[rowSums(umi.raw.subset) > 0,]
dim(umi.raw.subset_notzero)
dim(umi.raw.subset)

#get in correct formatting...
#construct sce so can use lognormcounts() function
umi.raw.subset.sce <- SingleCellExperiment(assays = list(counts = umi.raw.subset_notzero))

#add counts within each cell type for each donor ("pseudobulk") 
meta2$ctype <- as.character(meta2$ctype)
meta2$orig.ident <- as.character(meta2$orig.ident)
meta2$newclus <- paste0(meta2$orig.ident, ".",meta2$ctype)

Li_celltypes <- meta2$newclus[match(colnames(assays(umi.raw.subset.sce)$counts),rownames(meta2))]
umi.raw.subset.sce.bulked <- aggregateAcrossCells(umi.raw.subset.sce, ids=Li_celltypes)

save(umi.raw.subset.sce.bulked, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_bulked.rda")

#check summing
umitest <- assays(umi.raw.subset.sce)$counts
colnames(umitest) <- meta2$newclus[match(colnames(umitest),rownames(meta2))]
new2_umi2<- t(apply(umitest, 1, function(x) tapply(x, colnames(umitest), sum)))
bulked <-assays(umi.raw.subset.sce.bulked)$counts
all.equal(new2_umi2, bulked)
identical(new2_umi2, bulked, attrib.as.set=FALSE)
					
#calculate and divide bulked counts by scale factors
umi.raw.subset.sce.bulked.sf <- librarySizeFactors(umi.raw.subset.sce.bulked)
umi.raw.subset.sce.bulked.sfd <- logNormCounts(umi.raw.subset.sce.bulked, size_factors=umi.raw.subset.sce.bulked.sf, log=FALSE, center_size_factors = TRUE)

save(umi.raw.subset.sce.bulked.sfd, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_bulked_scalednorm.rda")

#log2 transform the bulked data in case its useful later, but we do not use this for calculating pSIs
umi.raw.subset.sce.bulked.sfd.log <- logNormCounts(umi.raw.subset.sce.bulked, size_factors=umi.raw.subset.sce.bulked.sf, log=TRUE, center_size_factors = TRUE)
save(umi.raw.subset.sce.bulked.sfd.log, file= "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_bulked_scalednormlog2.rda")

#for generating pSIs, get mean of biological reps (ie take mean of the three samples from the same cell type)
umi.raw.subset.sce.bulked.sfd.normcounts <- assays(umi.raw.subset.sce.bulked.sfd)$normcounts
Li_celltypes_abbr <- sub(".*\\.","",colnames(umi.raw.subset.sce.bulked.sfd.normcounts))

umi.means <- sumCountsAcrossCells(umi.raw.subset.sce.bulked.sfd.normcounts, ids = Li_celltypes_abbr, average = TRUE)

#check averaging
umitest2 <- umi.raw.subset.sce.bulked.sfd.normcounts
colnames(umitest2) <- Li_celltypes_abbr
new2<- t(apply(umitest2, 1, function(x) tapply(x, colnames(umitest2), mean)))
all.equal(new2, umi.means)
identical(new2, umi.means, attrib.as.set=FALSE)

#create filter data matrix where any values of 0 are replaced with NA
umi.filter <- umi.means
umi.filter[umi.filter == "0"] <- NA

#get pSIs
pSIs_umi_normcounts_pseudobulked_Li2018 <- specificity.index(umi.means, umi.filter, bts=50, p_max=1.0,e_min=0)

pSI.count(pSIs_umi_normcounts_pseudobulked_Li2018)


save(umi.means,umi.filter, pSIs_umi_normcounts_pseudobulked_Li2018, 
	file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda")

#################################################################################################################################

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)

#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")
gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Li.csv")		


#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")
gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
			
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Li.csv")		
				


##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#using Lake's frontal cortex data only
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers

#get data proper way
#process raw matrix
library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(SingleCellExperiment)
library(scater)
library(tibble)
library(dplyr)


setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/")

#################################################################################################################################

#get data for pSI calculation
fctx.raw.subset <- read.table("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")			
#note, Lake 2018's "raw" matrix file (above) is already subsetted to those nuclei that have passed QC

#drop any genes with no expression across all cells/nuclei
min(fctx.raw.subset)
sum(is.na(fctx.raw.subset))
totals <- rowSums(fctx.raw.subset)
min(totals)

fctx.raw.subset_notzero <- fctx.raw.subset[rowSums(fctx.raw.subset) > 0,]
dim(fctx.raw.subset_notzero)
dim(fctx.raw.subset)
#1833 genes were dropped

#get in correct formatting...

#construct sce so can use lognormcounts() function
fctx.raw.subset.sce <- SingleCellExperiment(assays = list(counts = fctx.raw.subset_notzero))

#using info from suppl table S1 from paper to match label to person and cell type
lake_meta2 <-read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/lake_meta2.csv",header=TRUE)
lake_meta2$OG <- as.character(paste0(lake_meta2$Identity, "_",lake_meta2$Names))
lake_meta2$Identityshort<-as.character(sub("[0-9].*","",lake_meta2$Identity))
lake_meta2$newIDs <- as.character(paste0(lake_meta2$Patient, ".", lake_meta2$Identityshort))
Lake_celltypes <- lake_meta2$newIDs[match(colnames(assays(fctx.raw.subset.sce)$counts),lake_meta2$OG)]

#add counts within each cell type for each donor ("pseudobulk") 
fctx.raw.subset.sce.bulked <- aggregateAcrossCells(fctx.raw.subset.sce, ids=Lake_celltypes)

save(fctx.raw.subset.sce.bulked, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_bulked.rda")

#check summing
umitest3 <- assays(fctx.raw.subset.sce)$counts
colnames(umitest3) <- lake_meta2$newIDs[match(colnames(umitest3),lake_meta2$OG)]
new3<- t(apply(umitest3, 1, function(x) tapply(x, colnames(umitest3), sum)))
bulked2 <-assays(fctx.raw.subset.sce.bulked)$counts
all.equal(new3, bulked2)
					
#calculate and divide bulked counts by scale factors
fctx.raw.subset.sce.bulked.sf <- librarySizeFactors(fctx.raw.subset.sce.bulked)
fctx.raw.subset.sce.bulked.sfd <- logNormCounts(fctx.raw.subset.sce.bulked, size_factors=fctx.raw.subset.sce.bulked.sf, log=FALSE, center_size_factors = TRUE)

save(fctx.raw.subset.sce.bulked.sfd, file="/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_bulked_scalednorm.rda")

#log2 transform the bulked data in case its useful later, but we do not use this for calculating pSIs
fctx.raw.subset.sce.bulked.sfd.log <- logNormCounts(fctx.raw.subset.sce.bulked, size_factors=fctx.raw.subset.sce.bulked.sf, log=TRUE, center_size_factors = TRUE)
save(fctx.raw.subset.sce.bulked.sfd.log, file= "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_bulked_scalednormlog2.rda")

#for generating pSIs, get mean of biological reps (ie take mean of the three samples from the same cell type)
fctx.raw.subset.sce.bulked.sfd.normcounts <- assays(fctx.raw.subset.sce.bulked.sfd)$normcounts
Lake_celltypes_abbr <- sub(".*\\.","",colnames(fctx.raw.subset.sce.bulked.sfd.normcounts))

fctx.means <- sumCountsAcrossCells(fctx.raw.subset.sce.bulked.sfd.normcounts, ids = Lake_celltypes_abbr, average = TRUE)

#check averaging
umitest4 <- fctx.raw.subset.sce.bulked.sfd.normcounts
colnames(umitest4) <- Lake_celltypes_abbr
new4<- t(apply(umitest4, 1, function(x) tapply(x, colnames(umitest4), mean)))
all.equal(new4, fctx.means)
identical(new4, fctx.means, attrib.as.set=FALSE)

#create filter data matrix where any values of 0 are replaced with NA
fctx.filter <- fctx.means
fctx.filter[fctx.filter == "0"] <- NA

#get pSIs
pSIs_fctx_normcounts_pseudobulked_Lake2018 <- specificity.index(fctx.means, fctx.filter, bts=50, p_max=1.0,e_min=0)

pSI.count(pSIs_fctx_normcounts_pseudobulked_Lake2018)

save(fctx.means,fctx.filter, pSIs_fctx_normcounts_pseudobulked_Lake2018, 
	file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda")

#################################################################################################################################

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda",verbose=TRUE)
					 
#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")
gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
			
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_amyg_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Lake.csv")		

#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")
gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_sacc_225_sym) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Lake.csv")		


###########################################################################################################################################################################
##eQTL genes: you can follow leo’s script for GO here: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/go_analysis.R to see what genes went in, but lets just 
	#start with the genes significantly enumerated on lines 105-106 (263/250/381)
############################################################################################################################################################################

#eQTL genes
#using gene symbols since that's the only input for CSEA

library('SummarizedExperiment')
library('data.table')
library('clusterProfiler')
library('GenomicRanges')
library('sessioninfo')
library('jaffelab')
library('pSI')
library(tibble)
library(dplyr)

setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr")

## Load the 881 eQTL results and get the unique gene IDs
eqtl_files <- dir(pattern = 'raggr_881_snps_.*')
names(eqtl_files) <- gsub('raggr_881_snps_|_eqtls_fdr01.csv', '', eqtl_files)

#subset into the different features
genes <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Gene",]
})

sapply(genes,dim)
#    amyg dlpfc sacc
#[1,] 2638  1886 4463
#[2,]   28    28   28


exons <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Exon",]
})

sapply(exons,dim)
#      amyg dlpfc  sacc
#[1,] 10370  6019 20014
#[2,]    28    28    28

junctions <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Jxn",]
})

sapply(junctions,dim)
#     amyg dlpfc sacc
#[1,] 6435  3927 8641
#[2,]   28    28   28


transcripts <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Tx",]
})

sapply(transcripts,dim)
#     amyg dlpfc sacc
#[1,] 3129  2837 5491
#[2,]   28    28   28


#get the unique, non-blank gene symbols for each feature

sig_genes <- lapply(genes,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_genes, length)
# amyg dlpfc  sacc 
#   84    89   144 

sig_exons <- lapply(exons,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_exons, length)
#amyg dlpfc  sacc 
#  132   126   208 

sig_junctions <- lapply(junctions,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_junctions, length)
#amyg dlpfc  sacc 
#  119   106   160 
sig_transcripts <- lapply(transcripts,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_transcripts, length)
#amyg dlpfc  sacc 
#   90   121   142 


#Subset the gene symbols by region

gene_amyg_eqtl <- sig_genes$amyg
exon_amyg_eqtl <- sig_exons$amyg
jxn_amyg_eqtl <- sig_junctions$amyg
tx_amyg_eqtl <- sig_transcripts$amyg

gene_dlpfc_eqtl <- sig_genes$dlpfc
exon_dlpfc_eqtl <- sig_exons$dlpfc
jxn_dlpfc_eqtl <- sig_junctions$dlpfc
tx_dlpfc_eqtl <- sig_transcripts$dlpfc

gene_sacc_eqtl <- sig_genes$sacc
exon_sacc_eqtl <- sig_exons$sacc
jxn_sacc_eqtl <- sig_junctions$sacc
tx_sacc_eqtl <- sig_transcripts$sacc


#Save objects as rdas and csvs
save(gene_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
save(exon_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
save(jxn_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
save(tx_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")

save(gene_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
save(exon_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
save(jxn_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
save(tx_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")

save(gene_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
save(exon_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
save(jxn_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
save(tx_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")

write.csv(gene_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_eQTL.csv")
write.csv(exon_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_Amyg_eQTL.csv")
write.csv(jxn_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_Amyg_eQTL.csv")
write.csv(tx_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_Amyg_eQTL.csv")

write.csv(gene_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_DLPFC_eQTL.csv")
write.csv(exon_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_DLPFC_eQTL.csv")
write.csv(jxn_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_DLPFC_eQTL.csv")
write.csv(tx_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_DLPFC_eQTL.csv")

write.csv(gene_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_eQTL.csv")
write.csv(exon_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_sACC_eQTL.csv")
write.csv(jxn_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_sACC_eQTL.csv")
write.csv(tx_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_sACC_eQTL.csv")



#all genes from all features together

## Load the 881 eQTL results and get the unique, not blank gene symbols from all features
eqtl_files <- dir(pattern = 'raggr_881_snps_.*')
names(eqtl_files) <- gsub('raggr_881_snps_|_eqtls_fdr01.csv', '', eqtl_files)

qtl <- lapply(eqtl_files, fread)

sig_genes_all <- lapply(qtl,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_genes_all, length)
# amyg dlpfc  sacc 
#  228   249   328 

gene_amyg_eqtl_all <- sig_genes_all$amyg
gene_dlpfc_eqtl_all <- sig_genes_all$dlpfc
gene_sacc_eqtl_all <- sig_genes_all$sacc

save(gene_amyg_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")
save(gene_dlpfc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")
save(gene_sacc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

write.csv(gene_amyg_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_eQTL_allFeatures.csv")
write.csv(gene_dlpfc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_DLPFC_eQTL_allFeatures.csv")
write.csv(gene_sacc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_eQTL_allFeatures.csv")


##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#CSEA using the pSI R package, built from source (http://genetics.wustl.edu/jdlab/psi_package/)
library(gdata)
library(pSI)
library(pSI.data)
data(mouse)
data(human)
library(tibble)
library(dplyr)

#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(mouse$psi.out,gene_amyg_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(mouse$psi.out,exon_amyg_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,jxn_amyg_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,tx_amyg_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,gene_amyg_eqtl_all,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(mouse$psi.out,gene_amyg_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,exon_amyg_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,jxn_amyg_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,tx_amyg_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Dougherty.csv")
write.csv(fisher.iteration(mouse$psi.out,gene_amyg_eqtl_all,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Dougherty.csv")				

#sACC, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(mouse$psi.out,gene_sacc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(mouse$psi.out,exon_sacc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,jxn_sacc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,tx_sacc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,gene_sacc_eqtl_all,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(mouse$psi.out,gene_sacc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,exon_sacc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,jxn_sacc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,tx_sacc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Dougherty.csv")
write.csv(fisher.iteration(mouse$psi.out,gene_sacc_eqtl_all,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Dougherty.csv")				

#DLPFC, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(mouse$psi.out,gene_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(mouse$psi.out,exon_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,jxn_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,tx_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,gene_dlpfc_eqtl_all,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(mouse$psi.out,gene_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,exon_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,jxn_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,tx_dlpfc_eqtl,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Dougherty.csv")
write.csv(fisher.iteration(mouse$psi.out,gene_dlpfc_eqtl_all,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Dougherty.csv")				

##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.

library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(tibble)
library(dplyr)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)

#Amyg, all features 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_amyg_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Li.csv")		

#DLPFC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_dlpfc_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_dlpfc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Li.csv")		

#sACC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,exon_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,jxn_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,tx_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,gene_sacc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Li.csv")		


##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#using Lake's frontal cortex data only
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers
library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(tibble)
library(dplyr)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda",verbose=TRUE)

#Amyg, all features 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_amyg_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_amyg_eqtl)	%>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_amyg_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Lake.csv")		

#DLPFC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_dlpfc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_dlpfc_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_dlpfc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Lake.csv")		

#sACC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda",verbose=TRUE)

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
				
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_sacc_eqtl) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
		
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_eqtl_all) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,exon_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,jxn_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,tx_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,gene_sacc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Lake.csv")		

#########################################################################################################################################################################################
##WGCNA modules: just which genes are in which modules to assign a cell type enrichment to each module: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna/analyze_results.R 
	#you can see how I got the genes for GO enrichment on lines 65-67
#########################################################################################################################################################################################
#this analysis is across brain regions (amygdala and sACC)
#this might be too many genes...for example, grey module has unique, non-blank 8817 genes
#other modules may have too few (~55)

library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
library(clusterProfiler)
library(readxl)
library(tibble)
library(dplyr)

setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/")

## load data
load("data/zandiHypde_bipolar_rseGene_n511.rda")
load("data/degradation_rse_BipSeq_BothRegions.rda")

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"),
                                levels = c("Control", "Bipolar"))

## add ancestry
load("genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

###########
# filter ##
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]
#number of gene symbols: 25136
#unique = 19201


#########################
## load wgcna output ####
#########################

## load
load("wgcna/rdas/constructed_network_signed_bicor.rda")

# get colors
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab,
        stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

#convert blank gene symbols to NA
length(rowData(rse_gene)$Symbol)
#25136
rowData(rse_gene)$Symbol[rowData(rse_gene)$Symbol == ""] <- NA
length(rowData(rse_gene)$Symbol)
#25136

#get unique, non-NA gene symbols
gList = split(rowData(rse_gene)$Symbol,
        factor(net$colorsLab,levels = colorDat$col))
sapply(gList, length)
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
sapply(gList, length)
gList = lapply(gList, function(x) as.character(unique(x[!is.na(x)])))
sapply(gList, length)
#the unique function only changed the grey module slightly, rest stayed same

#		 grey    turquoise         blue        brown       yellow        green 
#        8817         2147         2191         1489         1215         1043 
#         red        black         pink      magenta       purple  greenyellow 
#         518          273          216          217          196          180 
#         tan       salmon         cyan midnightblue    lightcyan       grey60 
#         137           93           75           75           77           59 
#  lightgreen  lightyellow    royalblue 
#          61           65           57 



#separate each module

genes_wgcna_grey <- gList$grey
genes_wgcna_turquoise <- gList$turquoise
genes_wgcna_blue <- gList$blue
genes_wgcna_brown <- gList$brown
genes_wgcna_yellow <- gList$yellow
genes_wgcna_green <- gList$green
genes_wgcna_red <- gList$red
genes_wgcna_black <- gList$black
genes_wgcna_pink <- gList$pink
genes_wgcna_magenta <- gList$magenta
genes_wgcna_purple <- gList$purple
genes_wgcna_greenyellow <- gList$greenyellow
genes_wgcna_tan <- gList$tan
genes_wgcna_salmon <- gList$salmon
genes_wgcna_cyan <- gList$cyan
genes_wgcna_midnightblue <- gList$midnightblue
genes_wgcna_lightcyan <- gList$lightcyan
genes_wgcna_grey60 <- gList$grey60
genes_wgcna_lightgreen <- gList$lightgreen
genes_wgcna_lightyellow <- gList$lightyellow
genes_wgcna_royalblue <- gList$royalblue

#spotcheck
length(genes_wgcna_grey)
#[1] 8817
length(genes_wgcna_turquoise)
#[1] 2147
length(genes_wgcna_pink)
#[1] 216
length(genes_wgcna_magenta)
#[1] 217
length(genes_wgcna_lightcyan)
#[1] 77
length(genes_wgcna_royalblue)
#[1] 57


#Save objects as rdas and csvs

save(genes_wgcna_grey, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
save(genes_wgcna_turquoise, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
save(genes_wgcna_blue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
save(genes_wgcna_brown, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
save(genes_wgcna_yellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
save(genes_wgcna_green, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
save(genes_wgcna_red, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
save(genes_wgcna_black, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
save(genes_wgcna_pink, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
save(genes_wgcna_magenta, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
save(genes_wgcna_purple, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
save(genes_wgcna_greenyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
save(genes_wgcna_tan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
save(genes_wgcna_salmon, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
save(genes_wgcna_cyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
save(genes_wgcna_midnightblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
save(genes_wgcna_lightcyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
save(genes_wgcna_grey60, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
save(genes_wgcna_lightgreen, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
save(genes_wgcna_lightyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
save(genes_wgcna_royalblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

write.csv(genes_wgcna_grey, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_grey.csv")
write.csv(genes_wgcna_turquoise, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_turquoise.csv")
write.csv(genes_wgcna_blue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_blue.csv")
write.csv(genes_wgcna_brown, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_brown.csv")
write.csv(genes_wgcna_yellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_yellow.csv")
write.csv(genes_wgcna_green, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_green.csv")
write.csv(genes_wgcna_red, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_red.csv")
write.csv(genes_wgcna_black, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_black.csv")
write.csv(genes_wgcna_pink, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_pink.csv")
write.csv(genes_wgcna_magenta, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_magenta.csv")
write.csv(genes_wgcna_purple, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_purple.csv")
write.csv(genes_wgcna_greenyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_greenyellow.csv")
write.csv(genes_wgcna_tan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_tan.csv")
write.csv(genes_wgcna_salmon, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_salmon.csv")
write.csv(genes_wgcna_cyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_cyan.csv")
write.csv(genes_wgcna_midnightblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_midnightblue.csv")
write.csv(genes_wgcna_lightcyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightcyan.csv")
write.csv(genes_wgcna_grey60, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_grey60.csv")
write.csv(genes_wgcna_lightgreen, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightgreen.csv")
write.csv(genes_wgcna_lightyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightyellow.csv")
write.csv(genes_wgcna_royalblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_royalblue.csv")


##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### CSEA

library(jaffelab)
library(SummarizedExperiment)
library(gdata)
library(pSI)
library(pSI.data)
data(mouse)
data(human)
library(tibble)
library(dplyr)

#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

#Calculation

fisher.iteration(mouse$psi.out,genes_wgcna_grey,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,genes_wgcna_turquoise,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_blue,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_brown,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_yellow,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_green,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_red,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_black,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_pink,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_magenta,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_purple,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_greenyellow,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_tan,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(mouse$psi.out,genes_wgcna_salmon,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_cyan,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_midnightblue,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_lightcyan,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_grey60,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_lightgreen,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_lightyellow,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(mouse$psi.out,genes_wgcna_royalblue,background="mouse.human",p.adjust=TRUE) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
#Save output/results

write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_grey,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Dougherty.csv")				
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_turquoise,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_blue,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_brown,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_yellow,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_green,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_red,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_black,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_pink,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_magenta,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_purple,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_greenyellow,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_tan,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_salmon,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_cyan,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_midnightblue,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_lightcyan,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_grey60,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_lightgreen,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_lightyellow,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Dougherty.csv")		
write.csv(fisher.iteration(mouse$psi.out,genes_wgcna_royalblue,background="mouse.human",p.adjust=TRUE), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Dougherty.csv")		

##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.

library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(tibble)
library(dplyr)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)

#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

#Calculation

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_grey) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_turquoise) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_blue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_brown) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_yellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_green) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_red) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_black) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_pink) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_magenta) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_purple) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_greenyellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_tan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_salmon) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_cyan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_midnightblue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightcyan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_grey60) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightgreen) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightyellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_royalblue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

#Save output/results

write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_grey), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Li.csv")				
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_turquoise), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_blue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_brown), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_yellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_green), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_red), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_black), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_pink), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_magenta), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_purple), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_greenyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_tan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_salmon), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_cyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_midnightblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightcyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_grey60), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightgreen), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_lightyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Li.csv")		
write.csv(fisher.iteration(pSIs_umi_normcounts_pseudobulked_Li2018,genes_wgcna_royalblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Li.csv")		


##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#using Lake's frontal cortex data only
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers

library(jaffelab)
library(SummarizedExperiment)
library(pSI)
library(tibble)
library(dplyr)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda",verbose=TRUE)

#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda",verbose=TRUE)

#Calculation

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_grey) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_turquoise) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_blue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_brown) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_yellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_green) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_red) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_black) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_pink)	%>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_magenta) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_purple) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_greenyellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_tan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_salmon) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_cyan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_midnightblue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightcyan) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_grey60) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightgreen) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightyellow) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")
	
fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_royalblue) %>%
	rownames_to_column(var = "feature") %>%
	filter_all(any_vars(. < 0.1)) %>%
	column_to_rownames(var = "feature")

#Save output/results

write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_grey), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Lake.csv")				
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_turquoise), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_blue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_brown), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_yellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_green), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_red), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_black), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_pink), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_magenta), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_purple), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_greenyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_tan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_salmon), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_cyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_midnightblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightcyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_grey60), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightgreen), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_lightyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Lake.csv")		
write.csv(fisher.iteration(pSIs_fctx_normcounts_pseudobulked_Lake2018,genes_wgcna_royalblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Lake.csv")		
	


