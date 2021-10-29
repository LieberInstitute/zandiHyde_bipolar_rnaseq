
library(variancePartition)
library(SingleCellExperiment)
library(tidyverse)
# library(limma)
# library(edgeR)
library(DeconvoBuddies)
library(jaffelab)
library(scuttle)
library(sessioninfo)
library(here)
library(scater)

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
marker_df <- read_csv(here("deconvolution", "bipseq-deconvolution-markers.csv"))
marker_genes <- marker_df$ensemblID

donor_region<- unique(colData(sce_pan)[, c("Sample","donor","region")])
table(donor_region$region)
# amy dlpfc   hpc   nac  sacc 
# 5     3     3     8     5 
## there are 24 samples

#### Pseudobulk sce samples ####
length(unique(sce_pan$sampleID))

pseuobulk_counts <- pseudobulk(sce_pan, cell_group_cols = c("sampleID", "cellType.Broad")) 
dim(pseuobulk_counts) # 23041 151
corner(pseuobulk_counts)

## pseudobulk pd
pseudobulk_pd <- as.data.frame(colData(sce_pan)) %>%
  select(sampleID, donor, region, protocol, cellType.Broad) %>%
  group_by(sampleID, donor, region, protocol, cellType.Broad) %>%
  summarise(n_cells = n())

rse_pseudobulk <- SummarizedExperiment(assays = list(counts = pseuobulk_counts),
                                       colData = DataFrame(pseudobulk_pd))

#### calcNormFactors approach ####
normFactors <- calcNormFactors(pseuobulk_counts)
pseuobulk_counts_cnf <- pseuobulk_counts*normFactors

## subset to marker genes
pseuobulk_counts_cnf <- pseuobulk_counts_cnf[marker_genes,]

message("Starting VarPart ", Sys.time())

# random effects as all are categorical
form <- ~(1|cellType.Broad) + (1|region) + (1|donor)

varPart <- fitExtractVarPartModel(exprObj = pseuobulk_counts_cnf, formula = form, data = pseudobulk_pd)
# save(varPart, file = "sce_refrence_variance_partition.rdata")

vp <- sortCols(varPart)
colMeans(vp)
# cellType.Broad         region          donor      Residuals 
#     0.45896885     0.03677204     0.01453426     0.48972486 

# cellType.Broad         region          donor      Residuals 
# 0.22778704     0.02136963     0.01949914     0.73134419 

vp_violin <- plotVarPart(vp) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(name = "Accent", n = 4)) + 
  labs(title = "Variance Partition for 225 Marker Genes - Normalized", 
       subtitle = "model: ~(1 | cellType.Broad) + (1 | region) + (1 | donor)")

ggsave(plot = vp_violin, filename = "plots/sc_refrence_vp_violin_norm.png", width = 10)

#### logNormCounts + getVarianceExplained approach ###

rse_pseudobulk <- logNormCounts(rse_pseudobulk)
vars <- getVarianceExplained(rse_pseudobulk, 
                             variables=c("cellType.Broad", "region", "donor"))
head(vars)
colMeans(vars)
# cellType.Broad         region          donor 
#      23.055392       3.450673       4.938777 

plot_var <- plotExplanatoryVariables(vars) + 
  theme_bw() +
  scale_x_continuous()
ggsave(plot_var, filename = "plots/sc_refrence_getVarExplain.png")

# sgejobs::job_single('sc_refrence_variance_partition', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript sc_refrence_variance_partition.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
