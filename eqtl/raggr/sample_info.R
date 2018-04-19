### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
pda = colData(rse_gene)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
pds = colData(rse_gene)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pdd = colData(rse_gene)

table(pds$PrimaryDx)
table(pda$PrimaryDx)
table(pdd$Dx)

table(pds$Sex)
table(pda$Sex)
table(pdd$Sex)

table(pds$PrimaryDx,pds$Sex)
table(pda$PrimaryDx,pda$Sex)
table(pdd$Dx,pdd$Sex)


library(VennDiagram)

venn.diagram(list(Amygdala = pda$BrNum, sACC = pds$BrNum, DLPFC = pdd$BrNum), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_sample_brains.png")

