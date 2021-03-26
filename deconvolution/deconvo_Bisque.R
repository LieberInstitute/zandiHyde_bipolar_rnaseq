
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(xbioc)
library(BisqueRNA)
library(tidyverse)
library(reshape2)
library(purrr)
library(compositions)
library(here)
library(sessioninfo)

#### Load Data ####
## Load rse_gene data
load(here("data", "zandiHypde_bipolar_rseGene_n513.rda"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/sce_filtered.Rdata", verbose = TRUE)

## marker data
load(here("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/marker_genes.Rdata"), verbose = TRUE)
all(marker_genes %in% rownames(rse_gene))
# [1] TRUE

#### create expression set ####
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene[marker_genes, ])$counts,
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce_pan[marker_genes,])$counts),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce_pan))[c("cellType.Broad", "cellType", "uniqueID","donor")]))

zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
exp_set_sce <- exp_set_sce[,zero_cell_filter]
  
#### Run Bisque ####
est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk, 
                                        sc.eset = exp_set_sce,
                                        cell.types = "cellType.Broad", 
                                        subject.names = "donor",
                                        use.overlap = FALSE)


est_prop$bulk.props <- t(est_prop$bulk.props)

est_prop$Est.prop.long <- melt(est_prop$bulk.props) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

est_prop$ilr <- ilr(est_prop$bulk.props)
colnames(est_prop$ilr) <- paste0("ilr_",1:ncol(est_prop$ilr))

est_prop_bisque <- est_prop

## Add long data and save
round(colMeans(est_prop_bisque$bulk.props),3)
# Astro Micro Oligo   OPC Excit Inhib 
# 0.105 0.089 0.570 0.081 0.077 0.078 

save(est_prop_bisque, file = here("deconvolution","est_prop_Bisque.Rdata"))

# sgejobs::job_single('deconvo_Bisque', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_Bisque.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
