library(SingleCellExperiment)
library(edgeR)
library(scran)
library(scuttle)

library(here)
library(sessioninfo)

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
table(sce_pan$region)
table(sce_pan$region, sce_pan$donor)

## filter to regions of intrest
sce_pan <- sce_pan[,sce_pan$region %in% c("amy","sacc")]
table(sce_pan$Sample)

summed <- aggregateAcrossCells(sce_pan, 
                               id=colData(sce_pan)[,c("cellType.Broad", "Sample")])
## Create DGRList
current <- summed[,summed$cellType.Broad == "Astro"]
y <- DGEList(counts(current), samples=colData(current))
dim(y)

## Romove samples w/ < 10 cells
discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded)

# Remove lowly expressed genes
keep <- filterByExpr(y, group=current$Sample)
# keep <- filterByExpr(y, group=current$region)
y <- y[keep,]
summary(keep)

## compute normalization
y <- calcNormFactors(y)
y$samples

## Check normalizations w/ plot
pdf("plots/sc_ref_norm.pdf")
par(mfrow=c(2,3))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}
dev.off()

## Model
design <- model.matrix(~factor(pool) + factor(tomato), y$samples)

#33 preform DE w/ pseudobulk data

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$celltype.mapped,
                            design=~factor(pool) + tomato,
                            coef="tomatoTRUE",
                            condition=summed.filt$tomato 
)

