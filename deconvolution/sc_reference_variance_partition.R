
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

donor_region<- unique(colData(sce_pan)[, c("Sample","donor","region")])
table(donor_region$region)
# amy dlpfc   hpc   nac  sacc 
# 5     3     3     8     5 
## there are 24 sampels

#### Pseudobulk sce samples ####
length(unique(sce_pan$sampleID))

pseuobulk_counts <- pseudobulk(sce_pan, cell_group_cols = c("sampleID", "cellType.Broad"))
dim(pseuobulk_counts)
corner(pseuobulk_counts)

## pseudobulk pd
pseudobulk_pd <- as.data.frame(colData(sce_pan)) %>%
  select(sampleID, donor, region, cellType.Broad) %>%
  group_by(sampleID, donor, region, cellType.Broad) %>%
  summarise(n_cells = n())

# mod <- model.matrix(~cellType.Broad, pd)
# gExpr<- calcNormFactors(sce_pan)
# vobjGenes <- voom(gExpr, mod)

message("Starting VarPart ", Sys.time())
form <- ~(1|cellType.Broad) + (1|region) + (1|donor)
# varPart <- fitExtractVarPartModel(exprObj = vobjGenes, formula = form, data = pd)
varPart <- fitExtractVarPartModel(exprObj = pseuobulk_counts, formula = form, data = pseudobulk_pd)
save(varPart, file = "sce_refrence_variance_partition.rdata")

vp <- sortCols(varPart)
vp_violin <- plotVarPart(vp) + labs(title = "mod")

ggsave(plot = vp_violin, filename = "plots/sc_refrence_vp_violin.png", width = 10)


# sgejobs::job_single('sc_refrence_variance_partition', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript sc_refrence_variance_partition.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
