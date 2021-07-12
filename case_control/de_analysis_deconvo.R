##### 

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)
library(purrr)
library(here)
# export _JAVA_OPTIONS="-Xms40g -Xmx60g"
library(xlsx)

## load data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load(here("data","degradation_rse_BipSeq_BothRegions.rda"), verbose = TRUE)

## deconvolution results
load(here("deconvolution","est_prop_Bisque.Rdata"),verbose = TRUE)

identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
identical(colnames(rse_gene), rownames(est_prop_bisque$bulk.props)) #TRUE

rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
				
## add ancestry 
load("../genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)


## filter
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

##############
## get qSVs ##
##############

modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 18
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 87%

# model w/o interaction to subset by region
modSep = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

modSep_prop <- cbind(modSep,est_prop_bisque$bulk.props[,c("Astro","Endo","Micro","Mural","Oligo","OPC","Tcell","Excit")])
modSep_ilr <- cbind(modSep,est_prop_bisque$ilr)

modSep_deconvo <- list(prop = modSep_prop, ilr = modSep_ilr)
#########################
## split back by region #
#########################

#### both ####
sACC_Index = which(colData(rse_gene)$BrainRegion == "sACC")
mod_sACC <- map(modSep_deconvo, ~cbind(.x[sACC_Index,], qSV_mat[sACC_Index, ]))

Amyg_Index = which(colData(rse_gene)$BrainRegion == "Amygdala")
mod_Amyg <- map(modSep_deconvo, ~cbind(.x[Amyg_Index,], qSV_mat[Amyg_Index, ]))

#### sACC ####
dge_sACC = DGEList(counts = assays(rse_gene[,sACC_Index])$counts, 
                   genes = rowData(rse_gene))
dge_sACC = calcNormFactors(dge_sACC)

outGene_sACC <- map(mod_sACC, function(m){

  vGene_sACC = voom(dge_sACC,m, plot=FALSE)
  fitGene_sACC = lmFit(vGene_sACC)
  eBGene_sACC = eBayes(fitGene_sACC)
  outGene_sACC = topTable(eBGene_sACC,coef=2,
                          p.value = 1,number=nrow(rse_gene))
  outGene_sACC = outGene_sACC[rownames(rse_gene),]
  return(outGene_sACC)
})

map_int(outGene_sACC, ~sum(.x$adj.P.Val < 0.05))
# prop  ilr 
# 318  379 

##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_gene[,Amyg_Index])$counts, 
	genes = rowData(rse_gene))
dge_Amyg = calcNormFactors(dge_Amyg)

outGene_Amyg <- map(mod_Amyg, function(m){

  vGene_Amyg = voom(dge_Amyg,m, plot=FALSE)
  
  fitGene_Amyg = lmFit(vGene_Amyg)
  eBGene_Amyg = eBayes(fitGene_Amyg)
  outGene_Amyg = topTable(eBGene_Amyg,coef=2,
                          p.value = 1,number=nrow(rse_gene))
  outGene_Amyg = outGene_Amyg[rownames(rse_gene),]
  return(outGene_Amyg)
})

map_int(outGene_Amyg, ~sum(.x$adj.P.Val < 0.05))

#### Export topTables to excel sheets ####
write.xlsx2(outGene_Amyg$prop, file = 'topTable_deconvolution.xlsx', sheetName='Amyg')
write.xlsx2(outGene_sACC$prop, file = 'topTable_deconvolution.xlsx', append = TRUE, sheetName='sACC')

write.xlsx2(outGene_Amyg$ilr, file = 'topTable_deconvolution_ilr.xlsx', sheetName='Amyg')
write.xlsx2(outGene_sACC$ilr, file = 'topTable_deconvolution_ilr.xlsx', append = TRUE, sheetName='sACC')

#### core output ####
nam = c("logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B", "deconvo_terms")

outGene_sACC_merged <- do.call("rbind",outGene_sACC)
outGene_sACC_merged$deconvo_terms <- ss(rownames(outGene_sACC_merged),"\\.")

outGene_Amyg_merged <- do.call("rbind",outGene_Amyg)
outGene_Amyg_merged$deconvo_terms <- ss(rownames(outGene_Amyg_merged),"\\.")

geneOut = cbind(outGene_Amyg_merged[,nam], outGene_sACC_merged[,nam])
colnames(geneOut) = paste0(colnames(geneOut), "_", rep(c("Amyg", "sACC"), each = length(nam)))
geneOut_deconvo <- geneOut

save(geneOut_deconvo, file = "bipolarControl_deStats_byRegion_qSVAjoint_deconvo.rda")

# sgejobs::job_single('de_analysis_deconvo', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript de_analysis_deconvo.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


# [1] "Reproducibility information:"
# [1] "2021-07-12 17:26:58 EDT"
# user  system elapsed 
# 289.244   4.990 301.307 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                 
# version  R version 4.0.4 RC (2021-02-08 r79975)
# os       CentOS Linux 7 (Core)                 
# system   x86_64, linux-gnu                     
# ui       X11                                   
# language (EN)                                  
# collate  en_US.UTF-8                           
# ctype    en_US.UTF-8                           
# tz       US/Eastern                            
# date     2021-07-12                            
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# annotate               1.68.0   2020-10-27 [2] Bioconductor                             
# AnnotationDbi          1.52.0   2020-10-27 [2] Bioconductor                             
# askpass                1.1      2019-01-13 [2] CRAN (R 4.0.3)                           
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.1    2020-12-09 [1] CRAN (R 4.0.4)                           
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocFileCache          1.14.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.1   2021-04-16 [2] Bioconductor                             
# BiocParallel         * 1.24.1   2020-11-06 [2] Bioconductor                             
# biomaRt                2.46.3   2021-02-09 [2] Bioconductor                             
# Biostrings             2.58.0   2020-10-27 [2] Bioconductor                             
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.0.3)                           
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.0.3)                           
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.0.4)                           
# blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.3)                           
# BSgenome               1.58.0   2020-10-27 [2] Bioconductor                             
# bumphunter             1.32.0   2020-10-27 [2] Bioconductor                             
# cachem                 1.0.4    2021-02-13 [2] CRAN (R 4.0.4)                           
# checkmate              2.0.0    2020-02-06 [2] CRAN (R 4.0.3)                           
# cli                    2.5.0    2021-04-26 [1] CRAN (R 4.0.4)                           
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.0.4)                           
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.0.4)                           
# colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)                           
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.0.3)                           
# curl                   4.3      2019-12-02 [2] CRAN (R 4.0.3)                           
# data.table             1.14.0   2021-02-21 [2] CRAN (R 4.0.4)                           
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.0.4)                           
# DelayedArray           0.16.3   2021-03-24 [2] Bioconductor                             
# derfinder              1.24.2   2020-12-18 [2] Bioconductor                             
# derfinderHelper        1.24.1   2020-12-18 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.0.3)                           
# downloader             0.4      2015-07-09 [2] CRAN (R 4.0.3)                           
# dplyr                  1.0.5    2021-03-05 [1] CRAN (R 4.0.4)                           
# edgeR                * 3.32.1   2021-01-14 [2] Bioconductor                             
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.0.4)                           
# fansi                  0.4.2    2021-01-15 [2] CRAN (R 4.0.3)                           
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.0.3)                           
# foreach                1.5.1    2020-10-15 [2] CRAN (R 4.0.3)                           
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.0.4)                           
# Formula                1.2-4    2020-10-16 [2] CRAN (R 4.0.3)                           
# genefilter           * 1.72.1   2021-01-21 [2] Bioconductor                             
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.7   2021-04-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicAlignments      1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicFeatures        1.42.3   2021-04-01 [2] Bioconductor                             
# GenomicFiles           1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# GEOquery               2.58.0   2020-10-27 [2] Bioconductor                             
# ggplot2                3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.0.3)                           
# Hmisc                  4.5-0    2021-02-28 [2] CRAN (R 4.0.4)                           
# hms                    1.0.0    2021-01-13 [2] CRAN (R 4.0.3)                           
# htmlTable              2.1.0    2020-09-16 [2] CRAN (R 4.0.3)                           
# htmltools              0.5.1.1  2021-01-22 [1] CRAN (R 4.0.4)                           
# htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.1   2020-12-12 [1] Bioconductor                             
# iterators              1.0.13   2020-10-15 [2] CRAN (R 4.0.3)                           
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 4.0.3)                           
# jsonlite               1.7.2    2020-12-09 [1] CRAN (R 4.0.3)                           
# knitr                  1.33     2021-04-24 [1] CRAN (R 4.0.4)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.4)                           
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.0.3)                           
# lifecycle              1.0.0    2021-02-15 [1] CRAN (R 4.0.4)                           
# limma                * 3.46.0   2020-10-27 [2] Bioconductor                             
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.3)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.4)                           
# MatrixGenerics       * 1.2.1    2021-01-30 [2] Bioconductor                             
# matrixStats          * 0.58.0   2021-01-29 [2] CRAN (R 4.0.3)                           
# memoise                2.0.0    2021-01-26 [2] CRAN (R 4.0.3)                           
# mgcv                 * 1.8-35   2021-04-18 [3] CRAN (R 4.0.4)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# nlme                 * 3.1-152  2021-02-04 [3] CRAN (R 4.0.4)                           
# nnet                   7.3-15   2021-01-24 [3] CRAN (R 4.0.4)                           
# openssl                1.4.4    2021-04-30 [1] CRAN (R 4.0.4)                           
# pillar                 1.6.0    2021-04-13 [1] CRAN (R 4.0.4)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.3)                           
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.3)                           
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.3)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# qvalue                 2.22.0   2020-10-27 [2] Bioconductor                             
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.0.3)                           
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.6    2021-01-15 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.0.4)                           
# readr                  1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# recount              * 1.16.1   2020-12-18 [2] Bioconductor                             
# rentrez                1.2.3    2020-11-10 [2] CRAN (R 4.0.3)                           
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rJava                  1.0-4    2021-04-29 [2] CRAN (R 4.0.4)                           
# rlang                  0.4.11   2021-04-30 [1] CRAN (R 4.0.4)                           
# rngtools               1.5      2020-01-23 [2] CRAN (R 4.0.3)                           
# rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.0.4)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# Rsamtools              2.6.0    2020-10-27 [2] Bioconductor                             
# RSQLite                2.2.7    2021-04-22 [2] CRAN (R 4.0.4)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rtracklayer            1.50.0   2020-10-27 [2] Bioconductor                             
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.0.4)                           
# sessioninfo            1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# survival               3.2-9    2021-03-14 [3] CRAN (R 4.0.4)                           
# sva                  * 3.38.0   2020-10-27 [2] Bioconductor                             
# tibble                 3.1.1    2021-04-18 [1] CRAN (R 4.0.4)                           
# tidyr                  1.1.3    2021-03-03 [2] CRAN (R 4.0.4)                           
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.0.4)                           
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.0.4)                           
# VariantAnnotation      1.36.0   2020-10-27 [2] Bioconductor                             
# vctrs                  0.3.8    2021-04-29 [1] CRAN (R 4.0.4)                           
# withr                  2.4.2    2021-04-18 [1] CRAN (R 4.0.4)                           
# xfun                   0.22     2021-03-11 [1] CRAN (R 4.0.4)                           
# xlsx                 * 0.6.5    2020-11-10 [2] CRAN (R 4.0.3)                           
# xlsxjars               0.6.1    2014-08-22 [2] CRAN (R 4.0.3)                           
# XML                    3.99-0.6 2021-03-16 [2] CRAN (R 4.0.4)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
# **** Job ends ****
#   Mon Jul 12 17:26:59 EDT 2021
