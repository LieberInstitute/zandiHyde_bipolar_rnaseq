
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
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)

n_genes <- list(`25` = 25, `50` = 50)

marker_stats <- map(n_genes, ~marker_stats %>% 
  filter(gene %in% statOut$ensemblID) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  dplyr::slice(1:.x))

map(marker_stats, ~.x%>% summarise(max(rank_ratio)))

tstats <- statOut[, grep("[f|t]_", colnames(statOut))]
fdrs <- statOut[, grep("adj.P.Val_", colnames(statOut))]
fdr_cut <- 0.05

marker_gene_list <- map(marker_stats,  function(ms){
  gene_list <- group_map(ms, ~pull(.x, gene))
  names(gene_list) <- levels(ms$cellType.target)
  return(gene_list)
})

source("gene_set_enrichment.R") 

gse <- map(marker_gene_list, ~gene_set_enrichment(gene_list = .x, modeling_results = statOut))
map(gse, ~.x %>% filter(Pval < 0.05))
# $`25`
#          OR         Pval        ID  test fdr_cut
# 1 47.263969 1.041614e-03 Amyg down Micro    0.05
# 2 22.204811 4.668394e-02 Amyg down Mural    0.05
# 3 22.204811 4.668394e-02 Amyg down Tcell    0.05
# 4 24.372437 9.364650e-08 sACC down Micro    0.05
# 5  8.465015 7.142635e-03 sACC down Mural    0.05
# 
# $`50`
#          OR         Pval        ID  test fdr_cut
# 1 19.293256 5.576193e-03   Amyg up  Endo    0.05
# 2 22.671134 4.125948e-03 Amyg down Micro    0.05
# 3  5.406839 8.245427e-03 sACC down  Endo    0.05
# 4 24.790286 3.073203e-14 sACC down Micro    0.05
# 5  5.406839 8.245427e-03 sACC down Mural    0.05

## save all values
walk2(gse, names(gse), ~write.csv(.x, file = here("case_control",paste0("gene_set_enrichment_cellType_markers",.y,".csv")), row.names = FALSE))

walk2(gse, names(gse), function(g, n){
  name <- here("case_control","plots",paste0("gene_set_enrichment_cellType_markers",n))
  t <- paste("OR: Top",n,"Cell Type Markers & FDR < 0.05")
  
  png(paste0(name,".png"))
  gene_set_enrichment_plot(g, ORcut = .2)
  title(t)
  dev.off()
  
  pdf(paste0(name, ".pdf"))
  gene_set_enrichment_plot(g, ORcut = .2)
  title(t)
  dev.off()
  
})


# sgejobs::job_single('de_cellType_marker_enrichment', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript de_cellType_marker_enrichment.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2021-07-12 17:37:37 EDT"
# > proc.time()
# user  system elapsed 
# 172.206  19.834 876.501 
# > options(width = 120)
# > session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.1.0 Patched (2021-05-18 r80330)
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
# package                * version  date       lib source                                   
# AnnotationDbi            1.54.1   2021-06-08 [2] Bioconductor                             
# AnnotationHub            3.0.1    2021-06-20 [2] Bioconductor                             
# assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)                           
# attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.0)                           
# backports                1.2.1    2020-12-09 [2] CRAN (R 4.1.0)                           
# beachmat                 2.8.0    2021-05-19 [2] Bioconductor                             
# beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.0)                           
# benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.0)                           
# benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.0)                           
# Biobase                * 2.52.0   2021-05-19 [2] Bioconductor                             
# BiocFileCache            2.0.0    2021-05-19 [2] Bioconductor                             
# BiocGenerics           * 0.38.0   2021-05-19 [2] Bioconductor                             
# BiocManager              1.30.16  2021-06-15 [2] CRAN (R 4.1.0)                           
# BiocNeighbors            1.10.0   2021-05-19 [1] Bioconductor                             
# BiocParallel             1.26.1   2021-07-04 [2] Bioconductor                             
# BiocSingular             1.8.1    2021-06-08 [1] Bioconductor                             
# BiocVersion              3.13.1   2021-03-19 [2] Bioconductor                             
# Biostrings               2.60.1   2021-06-06 [2] Bioconductor                             
# bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)                           
# bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)                           
# bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)                           
# blob                     1.2.1    2020-01-20 [2] CRAN (R 4.1.0)                           
# broom                    0.7.8    2021-06-24 [2] CRAN (R 4.1.0)                           
# bslib                    0.2.5.1  2021-05-18 [2] CRAN (R 4.1.0)                           
# cachem                   1.0.5    2021-05-15 [2] CRAN (R 4.1.0)                           
# cellranger               1.1.0    2016-07-27 [2] CRAN (R 4.1.0)                           
# cli                      3.0.0    2021-06-30 [2] CRAN (R 4.1.0)                           
# codetools                0.2-18   2020-11-04 [2] CRAN (R 4.1.0)                           
# colorout               * 1.2-2    2021-05-27 [1] Github (jalvesaq/colorout@79931fd)       
# colorspace               2.0-2    2021-06-24 [2] CRAN (R 4.1.0)                           
# config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.0)                           
# cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.0)                           
# crayon                   1.4.1    2021-02-08 [2] CRAN (R 4.1.0)                           
# curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)                           
# data.table               1.14.0   2021-02-21 [2] CRAN (R 4.1.0)                           
# DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)                           
# dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)                           
# DelayedArray             0.18.0   2021-05-19 [2] Bioconductor                             
# DelayedMatrixStats       1.14.0   2021-05-19 [2] Bioconductor                             
# desc                     1.3.0    2021-03-05 [2] CRAN (R 4.1.0)                           
# digest                   0.6.27   2020-10-24 [2] CRAN (R 4.1.0)                           
# dockerfiler              0.1.3    2019-03-19 [1] CRAN (R 4.1.0)                           
# doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)                           
# dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)                           
# dplyr                  * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)                           
# dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.0)                           
# DropletUtils             1.12.1   2021-06-01 [1] Bioconductor                             
# DT                       0.18     2021-04-14 [2] CRAN (R 4.1.0)                           
# edgeR                    3.34.0   2021-05-19 [2] Bioconductor                             
# ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)                           
# ExperimentHub            2.0.0    2021-05-19 [2] Bioconductor                             
# fansi                    0.5.0    2021-05-25 [2] CRAN (R 4.1.0)                           
# fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)                           
# fields                   12.5     2021-06-25 [2] CRAN (R 4.1.0)                           
# filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)                           
# forcats                * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)                           
# foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)                           
# fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)                           
# gargle                   1.2.0    2021-07-02 [2] CRAN (R 4.1.0)                           
# generics                 0.1.0    2020-10-31 [2] CRAN (R 4.1.0)                           
# GenomeInfoDb           * 1.28.1   2021-07-01 [2] Bioconductor                             
# GenomeInfoDbData         1.2.6    2021-05-11 [2] Bioconductor                             
# GenomicRanges          * 1.44.0   2021-05-19 [2] Bioconductor                             
# ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.0)                           
# ggplot2                * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)                           
# glue                     1.4.2    2020-08-27 [2] CRAN (R 4.1.0)                           
# golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.0)                           
# googledrive              2.0.0    2021-07-08 [2] CRAN (R 4.1.0)                           
# gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)                           
# gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)                           
# haven                    2.4.1    2021-04-23 [2] CRAN (R 4.1.0)                           
# HDF5Array                1.20.0   2021-05-19 [2] Bioconductor                             
# here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)                           
# hms                      1.1.0    2021-05-17 [2] CRAN (R 4.1.0)                           
# htmltools                0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)                           
# htmlwidgets              1.5.3    2020-12-10 [2] CRAN (R 4.1.0)                           
# httpuv                   1.6.1    2021-05-07 [2] CRAN (R 4.1.0)                           
# httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)                           
# interactiveDisplayBase   1.30.0   2021-05-19 [2] Bioconductor                             
# IRanges                * 2.26.0   2021-05-19 [2] Bioconductor                             
# irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)                           
# iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)                           
# jaffelab               * 0.99.31  2021-05-27 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)                           
# jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)                           
# KEGGREST                 1.32.0   2021-05-19 [2] Bioconductor                             
# knitr                    1.33     2021-04-24 [2] CRAN (R 4.1.0)                           
# later                    1.2.0    2021-04-23 [2] CRAN (R 4.1.0)                           
# lattice                  0.20-44  2021-05-02 [3] CRAN (R 4.1.0)                           
# lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)                           
# lifecycle                1.0.0    2021-02-15 [2] CRAN (R 4.1.0)                           
# limma                    3.48.1   2021-06-24 [2] Bioconductor                             
# locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)                           
# lubridate                1.7.10   2021-02-26 [2] CRAN (R 4.1.0)                           
# magick                   2.7.2    2021-05-02 [2] CRAN (R 4.1.0)                           
# magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)                           
# maps                     3.3.0    2018-04-03 [2] CRAN (R 4.1.0)                           
# Matrix                   1.3-4    2021-06-01 [3] CRAN (R 4.1.0)                           
# MatrixGenerics         * 1.4.0    2021-05-19 [2] Bioconductor                             
# matrixStats            * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)                           
# memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)                           
# mime                     0.11     2021-06-23 [2] CRAN (R 4.1.0)                           
# modelr                   0.1.8    2020-05-19 [2] CRAN (R 4.1.0)                           
# munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)                           
# pillar                   1.6.1    2021-05-16 [2] CRAN (R 4.1.0)                           
# pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)                           
# pkgload                  1.2.1    2021-04-06 [2] CRAN (R 4.1.0)                           
# plotly                   4.9.4.1  2021-06-18 [2] CRAN (R 4.1.0)                           
# png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)                           
# Polychrome               1.2.6    2020-11-11 [1] CRAN (R 4.1.0)                           
# promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)                           
# purrr                  * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)                           
# R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)                           
# R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)                           
# R.utils                  2.10.1   2020-08-26 [2] CRAN (R 4.1.0)                           
# R6                       2.5.0    2020-10-28 [2] CRAN (R 4.1.0)                           
# rafalib                * 1.0.0    2015-08-09 [1] CRAN (R 4.1.0)                           
# rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)                           
# RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)                           
# Rcpp                     1.0.7    2021-07-07 [2] CRAN (R 4.1.0)                           
# RCurl                    1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)                           
# readr                  * 1.4.0    2020-10-05 [2] CRAN (R 4.1.0)                           
# readxl                   1.3.1    2019-03-13 [2] CRAN (R 4.1.0)                           
# remotes                  2.4.0    2021-06-02 [2] CRAN (R 4.1.0)                           
# reprex                   2.0.0    2021-04-02 [2] CRAN (R 4.1.0)                           
# rhdf5                    2.36.0   2021-05-19 [2] Bioconductor                             
# rhdf5filters             1.4.0    2021-05-19 [2] Bioconductor                             
# Rhdf5lib                 1.14.2   2021-07-06 [2] Bioconductor                             
# rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)                           
# rlang                    0.4.11   2021-04-30 [2] CRAN (R 4.1.0)                           
# roxygen2                 7.1.1    2020-06-27 [2] CRAN (R 4.1.0)                           
# rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)                           
# RSQLite                  2.2.7    2021-04-22 [2] CRAN (R 4.1.0)                           
# rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)                           
# rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.0)                           
# rvest                    1.0.0    2021-03-09 [2] CRAN (R 4.1.0)                           
# S4Vectors              * 0.30.0   2021-05-19 [2] Bioconductor                             
# sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)                           
# ScaledMatrix             1.0.0    2021-05-19 [1] Bioconductor                             
# scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)                           
# scater                   1.20.1   2021-06-15 [1] Bioconductor                             
# scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.0)                           
# scuttle                  1.2.0    2021-05-19 [1] Bioconductor                             
# segmented                1.3-4    2021-04-22 [1] CRAN (R 4.1.0)                           
# sessioninfo            * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)                           
# shiny                    1.6.0    2021-01-25 [2] CRAN (R 4.1.0)                           
# shinyWidgets             0.6.0    2021-03-15 [1] CRAN (R 4.1.0)                           
# SingleCellExperiment   * 1.14.1   2021-05-21 [2] Bioconductor                             
# spam                     2.7-0    2021-06-25 [2] CRAN (R 4.1.0)                           
# sparseMatrixStats        1.4.0    2021-05-19 [2] Bioconductor                             
# SpatialExperiment      * 1.2.1    2021-06-10 [1] Bioconductor                             
# spatialLIBD            * 1.4.0    2021-05-20 [1] Bioconductor                             
# stringi                  1.6.2    2021-05-17 [2] CRAN (R 4.1.0)                           
# stringr                * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)                           
# SummarizedExperiment   * 1.22.0   2021-05-19 [2] Bioconductor                             
# testthat                 3.0.4    2021-07-01 [2] CRAN (R 4.1.0)                           
# tibble                 * 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)                           
# tidyr                  * 1.1.3    2021-03-03 [2] CRAN (R 4.1.0)                           
# tidyselect               1.1.1    2021-04-30 [2] CRAN (R 4.1.0)                           
# tidyverse              * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)                           
# usethis                  2.0.1    2021-02-10 [2] CRAN (R 4.1.0)                           
# utf8                     1.2.1    2021-03-12 [2] CRAN (R 4.1.0)                           
# vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)                           
# vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.0)                           
# viridis                  0.6.1    2021-05-11 [2] CRAN (R 4.1.0)                           
# viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)                           
# withr                    2.4.2    2021-04-18 [2] CRAN (R 4.1.0)                           
# xfun                     0.24     2021-06-15 [2] CRAN (R 4.1.0)                           
# xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)                           
# xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)                           
# XVector                  0.32.0   2021-05-19 [2] Bioconductor                             
# yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)                           
# zlibbioc                 1.38.0   2021-05-19 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
