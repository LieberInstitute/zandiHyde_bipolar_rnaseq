
library(variancePartition)
library(SingleCellExperiment)
library(tidyverse)
library(limma)
library(edgeR)
library(DeconvoBuddies)
library(jaffelab)
library(sessioninfo)
library(here)

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)
marker_stats_filter <- marker_stats %>%
  filter(gene %in% rownames(rse_gene)) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  slice(1:25)

marker_genes <- marker_stats_filter$gene

donor_region<- unique(colData(sce_pan)[, c("Sample","donor","region")])
table(donor_region$region)
# amy dlpfc   hpc   nac  sacc 
# 5     3     3     8     5 
## there are 24 samples

#### Pseudobulk sce samples ####
length(unique(sce_pan$sampleID))

pseuobulk_counts <- pseudobulk(sce_pan[marker_genes,], cell_group_cols = c("sampleID", "cellType.Broad"))
dim(pseuobulk_counts) # 225 151
corner(pseuobulk_counts)

## pseudobulk pd
pseudobulk_pd <- as.data.frame(colData(sce_pan)) %>%
  select(sampleID, donor, region, protocol, cellType.Broad) %>%
  group_by(sampleID, donor, region, protocol, cellType.Broad) %>%
  summarise(n_cells = n())

rse_pseudobulk <- SummarizedExperiment(assays = list(counts = pseuobulk_counts),
                                       colData = DataFrame(pseudobulk_pd))


message("Starting VarPart ", Sys.time())
# linear mixed model
form <- ~(1|cellType.Broad) + (1|region) + (1|donor)
varPart <- fitExtractVarPartModel(exprObj = pseuobulk_counts, formula = form, data = pseudobulk_pd)
save(varPart, file = "sce_refrence_variance_partition.rdata")

vp <- sortCols(varPart)
colMeans(vp)
# cellType.Broad         region          donor      Residuals 
#     0.45896885     0.03677204     0.01453426     0.48972486 

vp_violin <- plotVarPart(vp) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(name = "Accent", n = 4)) + 
  labs(title = "Variance Partition for 225 Marker Genes", 
       subtitle = "model: ~(1 | cellType.Broad) + (1 | region) + (1 | donor)")

ggsave(plot = vp_violin, filename = "plots/sc_refrence_vp_violin.png", width = 10)


# sgejobs::job_single('sc_refrence_variance_partition', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript sc_refrence_variance_partition.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
