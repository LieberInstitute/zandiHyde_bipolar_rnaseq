
library(variancePartition)
library(SummarizedExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(jaffelab)
library(edgeR)
library(limma)
library(voom)
library(sessioninfo)
library(here)
# library(scater)

## load data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)

rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
                     levels = c("Control", "Bipolar"))

marker_df <- read_csv(here("deconvolution", "bipseq-deconvolution-markers.csv"))
marker_genes <- marker_df$gencodeID


gExpr <- DGEList(assays(rse_gene)$counts)
gExpr <- calcNormFactors(gExpr)
vobjGenes <- voom(gExpr)

vobjGenes <- vobjGenes[marker_genes,] ## subset to 225 marker genes
dim(vobjGenes) 
# [1] 225 511

## Load other data needed for models
load(here("case_control","qSV_mat.rda"), verbose = TRUE)
load(here("deconvolution","est_prop_Bisque.Rdata"), verbose = TRUE)
load(here("genotype_data","zandiHyde_bipolar_MDS_n511.rda"), verbose = TRUE)
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)

info <- as.data.frame(cbind(colData(rse_gene),mds, qSV_mat, est_prop_bisque$bulk.props) )

## 1. Simple Model based on what we did for SC refrence
## form <- ~(1|cellType.Broad) + (1|region) + (1|donor)
formSimple <- ~ (1|Dx) + (1|BrainRegion) + (1|BrNum) + Astro + Endo + Micro + 
  Mural + Oligo + OPC + Tcell + Excit

varPart_simple <- fitExtractVarPartModel( vobjGenes, formSimple, info)
cellFractions <- rowSums(varPart_simple[c("Astro", "Endo", "Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit")])
varPart_simple_sum <- cbind(cellFractions, varPart_simple[c("Dx", "BrainRegion", "BrNum", "Residuals")]) 
varPart_simple_sum <- sortCols(varPart_simple_sum)

## 2. Final model we ran on bulk data, based on modJoint
formJoint <- ~Dx + BrainRegion + Dx:BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
  PC1 + PC2 + PC3 + PC4 +  PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + 
  PC13 + PC14 + PC15 + PC16 + PC17 + PC18

varPart_joint <- fitExtractVarPartModel( vobjGenes, formJoint, info)
qSVs <- rowSums(varPart_joint[,grepl("^PC", colnames(varPart_joint))])
varPart_joint_sum <- cbind(qSVs, varPart_joint[,!grepl("^PC", colnames(varPart_joint))])
varPart_joint_sum <- sortCols(varPart_joint_sum)
colMeans(varPart_joint_sum)

save(varPart_simple, varPart_simple_sum, varPart_joint, varPart_joint_sum,file = "bulk_varPart.Rdata")

#### Plot ####
violin_simple <- plotVarPart(varPart_simple_sum, col = rep("white", ncol(varPart_simple_sum))) +
  labs(y = "Variance explained (%)") +
  theme(text = element_text(size=15))
ggsave(violin_simple, filename = "plots/bulk_vp_violin_simple.png", width = 4)

violin_joint <- plotVarPart(varPart_joint_sum, col = rep("white", ncol(varPart_joint_sum)))+
  labs(y = "Variance explained (%)") +
  theme(text = element_text(size=15)) 
ggsave(violin_joint, filename = "plots/bulk_vp_violin_joint.png", width = 10, height = 7.5) 

# sgejobs::job_single('bulk_variance_partition', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript bulk_variance_partition.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
