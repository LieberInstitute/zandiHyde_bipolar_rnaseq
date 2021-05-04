
library(tidyverse)
library(spatialLIBD)
library(SummarizedExperiment)
library(jaffelab)
library(here)
library(sessioninfo)

## Load Data
# load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

statOut <- statOut %>% 
  as.data.frame() %>% 
  filter(Type == "Gene") %>%
  mutate(ensemblID = ss(gencodeID,"\\."))

dim(statOut)

## marker data
load(here("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/marker_stats.Rdata"), verbose = TRUE)

marker_stats_50 <- marker_stats %>% 
  group_by(cellType.target) %>%
  filter(rank_ratio <= 50)

all(marker_stats_50$gene %in% statOut$ensemblID)
# [1] TRUE

tstats <- statOut[, grep("[f|t]_", colnames(statOut))]
fdrs <- statOut[, grep("adj.P.Val_", colnames(statOut))]
fdr_cut <- 0.05

marker_gene_list <- marker_stats_50 %>% 
  group_map(~pull(.x, gene))
names(marker_gene_list) <- levels(marker_stats_50$cellType.target)


source("gene_set_enrichment.R") 

gse <- gene_set_enrichment(gene_list = marker_gene_list, modeling_results = statOut)
gse %>% filter(Pval < 0.05)
#         OR         Pval   test Group direction fdr_cut
# 1 24.79029 3.073203e-14 t_sACC Micro      down    0.05

gse %>% filter(Pval < .99) %>% summarize(max(Pval))

gse$ID <- gsub("_"," ",gsub("t_","",gse$ID))

png(here("case_control","plots","gene_set_enrichment.png"))
gene_set_enrichment_plot(gse, ORcut = .2)
title("OR: Top 50 Cell Type Markers & FDR < 0.05")
dev.off()

pdf(here("case_control","plots","gene_set_enrichment.pdf"))
gene_set_enrichment_plot(gse, ORcut = .2)
title("OR: Top 50 Cell Type Markers & FDR < 0.05")
dev.off()

# sgejobs::job_single('de_cellType_marker_enrichment', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript de_cellType_marker_enrichment.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
