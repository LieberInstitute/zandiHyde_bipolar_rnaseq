
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
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
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

# [1] "Reproducibility information:"
# [1] "2021-03-26 16:20:17 EDT"
# user  system elapsed 
# 40.506   2.638  44.138 
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
# date     2021-03-26                            
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date       lib source                                   
# AnnotationDbi        * 1.52.0     2020-10-27 [2] Bioconductor                             
# assertthat             0.2.1      2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.1      2020-12-09 [1] CRAN (R 4.0.4)                           
# bayesm                 3.1-4      2019-10-15 [1] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0     2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0     2020-10-27 [2] Bioconductor                             
# BiocManager            1.30.10    2019-11-16 [2] CRAN (R 4.0.3)                           
# BisqueRNA            * 1.0.4      2020-07-18 [1] CRAN (R 4.0.3)                           
# bit                    4.0.4      2020-08-04 [2] CRAN (R 4.0.3)                           
# bit64                  4.0.5      2020-08-30 [2] CRAN (R 4.0.3)                           
# bitops                 1.0-6      2013-08-17 [2] CRAN (R 4.0.3)                           
# blob                   1.2.1      2020-01-20 [2] CRAN (R 4.0.3)                           
# broom                  0.7.5      2021-02-19 [2] CRAN (R 4.0.4)                           
# cachem                 1.0.4      2021-02-13 [2] CRAN (R 4.0.4)                           
# cellranger             1.1.0      2016-07-27 [2] CRAN (R 4.0.3)                           
# checkmate              2.0.0      2020-02-06 [2] CRAN (R 4.0.3)                           
# cli                    2.3.1      2021-02-23 [1] CRAN (R 4.0.4)                           
# codetools              0.2-18     2020-11-04 [3] CRAN (R 4.0.4)                           
# colorspace             2.0-0      2020-11-11 [2] CRAN (R 4.0.3)                           
# compositions         * 2.0-1      2021-01-08 [1] CRAN (R 4.0.3)                           
# crayon                 1.4.1      2021-02-08 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.1      2021-01-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.1.0      2021-02-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.3     2021-03-24 [2] Bioconductor                             
# DEoptimR               1.0-8      2016-11-19 [2] CRAN (R 4.0.3)                           
# digest                 0.6.27     2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                * 1.0.5      2021-03-05 [1] CRAN (R 4.0.4)                           
# ellipsis               0.3.1      2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.2      2021-01-15 [2] CRAN (R 4.0.3)                           
# fastmap                1.1.0      2021-01-25 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.1      2021-01-27 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0      2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0      2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.4     2021-03-10 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4      2020-11-30 [2] Bioconductor                             
# GenomicRanges        * 1.42.0     2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3      2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2      2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1      2020-05-05 [1] CRAN (R 4.0.3)                           
# gtable                 0.3.0      2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1      2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.0.3)                           
# hms                    1.0.0      2021-01-13 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2      2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.1     2020-12-12 [1] Bioconductor                             
# jaffelab             * 0.99.30    2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jsonlite               1.7.2      2020-12-09 [1] CRAN (R 4.0.3)                           
# lattice                0.20-41    2020-04-02 [3] CRAN (R 4.0.4)                           
# lifecycle              1.0.0      2021-02-15 [1] CRAN (R 4.0.4)                           
# limma                  3.46.0     2020-10-27 [2] Bioconductor                             
# limSolve               1.5.6      2019-11-14 [1] CRAN (R 4.0.3)                           
# lpSolve                5.6.15     2020-01-24 [1] CRAN (R 4.0.3)                           
# lubridate              1.7.10     2021-02-26 [1] CRAN (R 4.0.4)                           
# magrittr               2.0.1      2020-11-17 [2] CRAN (R 4.0.3)                           
# MASS                   7.3-53.1   2021-02-12 [3] CRAN (R 4.0.4)                           
# Matrix                 1.3-2      2021-01-06 [3] CRAN (R 4.0.4)                           
# MatrixGenerics       * 1.2.1      2021-01-30 [2] Bioconductor                             
# matrixStats          * 0.58.0     2021-01-29 [2] CRAN (R 4.0.3)                           
# memoise                2.0.0      2021-01-26 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8      2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0      2018-06-12 [2] CRAN (R 4.0.3)                           
# pillar                 1.5.1      2021-03-05 [1] CRAN (R 4.0.4)                           
# pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.0.3)                           
# pkgmaker               0.32.2.900 2020-11-02 [1] Github (renozao/pkgmaker@51b7207)        
# plyr                   1.8.6      2020-03-03 [2] CRAN (R 4.0.3)                           
# ps                     1.6.0      2021-02-28 [1] CRAN (R 4.0.4)                           
# purrr                * 0.3.4      2020-04-17 [1] CRAN (R 4.0.3)                           
# quadprog               1.5-8      2019-11-20 [2] CRAN (R 4.0.3)                           
# R6                     2.5.0      2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.0.3)                           
# RColorBrewer           1.1-2      2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.6      2021-01-15 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.3   2021-03-16 [2] CRAN (R 4.0.4)                           
# readr                * 1.4.0      2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1      2019-03-13 [2] CRAN (R 4.0.3)                           
# registry               0.5-1      2019-03-05 [2] CRAN (R 4.0.3)                           
# reprex                 1.0.0      2021-01-27 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4      2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                  0.4.10     2020-12-30 [1] CRAN (R 4.0.4)                           
# robustbase             0.93-7     2021-01-04 [2] CRAN (R 4.0.3)                           
# rprojroot              2.0.2      2020-11-15 [2] CRAN (R 4.0.3)                           
# RSQLite                2.2.4      2021-03-12 [2] CRAN (R 4.0.4)                           
# rstudioapi             0.13       2020-11-12 [2] CRAN (R 4.0.3)                           
# rvest                  1.0.0      2021-03-09 [2] CRAN (R 4.0.4)                           
# S4Vectors            * 0.28.1     2020-12-09 [2] Bioconductor                             
# scales                 1.1.1      2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-3      2021-03-08 [1] CRAN (R 4.0.4)                           
# sessioninfo          * 1.1.1      2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0     2020-10-27 [2] Bioconductor                             
# stringi                1.5.3      2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0      2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0     2020-10-27 [2] Bioconductor                             
# tensorA                0.36.2     2020-11-19 [1] CRAN (R 4.0.3)                           
# tibble               * 3.1.0      2021-02-25 [1] CRAN (R 4.0.4)                           
# tidyr                * 1.1.3      2021-03-03 [2] CRAN (R 4.0.4)                           
# tidyselect             1.1.0      2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0      2019-11-21 [1] CRAN (R 4.0.3)                           
# utf8                   1.2.1      2021-03-12 [2] CRAN (R 4.0.4)                           
# vctrs                  0.3.6      2020-12-17 [1] CRAN (R 4.0.4)                           
# withr                  2.4.1      2021-01-26 [1] CRAN (R 4.0.3)                           
# xbioc                * 0.1.19     2020-11-02 [1] Github (renozao/xbioc@1354168)           
# xml2                   1.3.2      2020-04-23 [2] CRAN (R 4.0.3)                           
# xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0     2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0     2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
# **** Job ends ****
#   Fri Mar 26 16:20:18 EDT 2021
# 
