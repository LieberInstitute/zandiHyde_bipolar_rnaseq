library(SingleCellExperiment)
library(edgeR)
library(scran)
library(scuttle)
library(here)
library(sessioninfo)
library(purrr)

## Use V1 sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)


## filter to regions of intrest
table(sce_pan$region)
table(sce_pan$region, sce_pan$donor)
sce_pan <- sce_pan[,sce_pan$region %in% c("amy","sacc")]

table(sce_pan$sampleID, sce_pan$cellType.Broad)

summed <- aggregateAcrossCells(sce_pan, 
                               id=colData(sce_pan)[,c("cellType.Broad", "sampleID")])

colData(summed)

#33 preform DE w/ pseudobulk data
summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$cellType.Broad,
                            design= ~region + donor,
                            coef="regionsacc",
                            condition=summed.filt$region
)

map(list(de.results),~sum(.x$FDR < 0.05))
sum(de.results$Astro$FDR < 0.05, na.rm = TRUE)
sum(de.results$Oligo$FDR < 0.05, na.rm = TRUE)
