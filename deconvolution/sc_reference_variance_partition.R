
library(variancePartition)
library(SingleCellExperiment)
library(purrr)
library(limma)
library(edgeR)
library(sessioninfo)
library(here)

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)

donor_region<- unique(colData(sce_pan)[, c("Sample","donor","region")])
table(donor_region$region)
# amy dlpfc   hpc   nac  sacc 
# 5     3     3     8     5 

pd <- colData(sce_pan)

mod <- model.matrix(~cellType.Broad, pd)
gExpr<- calcNormFactors(sce_pan)
vobjGenes <- voom(gExpr, mod)

message("VarPart Amyg")
form <- ~(1|cellType.Broad) + (1|region) + (1|donor)
varPart <- fitExtractVarPartModel(exprObj = vobjGenes, formula = form, data = pd)
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
