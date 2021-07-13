
library(SummarizedExperiment)
library(jaffelab)
library(matrixStats)
library(recount)
library(here)
library(sessioninfo)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)

## load data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load(here("data","zandiHypde_bipolar_rseExon_n511.rda"), verbose = TRUE)
load(here("data","zandiHypde_bipolar_rseJxn_n511.rda"), verbose = TRUE)
load(here("data","zandiHypde_bipolar_rseTx_n511.rda"), verbose = TRUE)

rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"),
                     levels = c("Control", "Bipolar"))

assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25
rse_gene = rse_gene[geneIndex,]

dim(rse_gene)
# [1] 25136   511

## add ancestry
load(here("genotype_data","zandiHyde_bipolar_MDS_n511.rda"), verbose = TRUE)
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## modJoint
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                        data=colData(rse_gene))

modJoint <- modJoint[,c(1:3,14,4:13)]
colnames(modJoint)

assays(rse_exon)$rpkm <- getRPKM(rse_exon)
assays(rse_jxn)$rpm <- getRPM(rse_jxn)

data <- list(Gene = rse_gene, Exon = rse_exon, Jxn = rse_jxn, Tx = rse_tx)
##Run cleaningY
# cleaningY(log2(RPKM + 1)) keeping Dx and Region effects and removing everything else (except the intercept)
expr_gene <- cleaningY(log2(assays(rse_gene)$rpkm + 1), mod = modJoint, P = 4)
expr_exon <- cleaningY(log2(assays(rse_exon)$rpkm + 1), mod = modJoint, P = 4)
expr_jxn <- cleaningY(log2(assays(rse_jxn)$rpm + 1), mod = modJoint, P = 4)
expr_tx <- cleaningY(log2(assays(rse_tx)$tpm + 1), mod = modJoint, P = 4)

expr <- list(Gene = expr_gene, Exon = expr_exon, Jxn = expr_jxn, Tx = expr_tx)

#### Box Plots ####
region_colors <- list(Amygdala = "#FFFF1F",
                      sACC = "#8EB0F6")

# it’s going to be 4 boxes (region * Dx), with the rowSds() output from the

# regionXdx <- cross2(levels(rse_gene$BrainRegion),levels(rse_gene$Dx))
# names(regionXdx) <- map(regionXdx, ~paste(.x[[1]],.x[[2]]))
# 
# table(rse_gene$BrainRegion, rse_gene$Dx)
# 
# gene_sd <- map_dfc(regionXdx, function(rXd){
#   expr_gene_subset <- expr_gene[,rse_gene$BrainRegion == rXd[[1]] & rse_gene$Dx == rXd[[2]]]
#   return(rowSds(expr_gene_subset))
# })
# 
# gene_sd_long <- gene_sd %>%
#   add_column(gene = rownames(expr_gene)) %>%
#   pivot_longer(!gene, names_to = "regionXdx", values_to = "gene_sd") %>%
#   mutate(regionXdx2 = regionXdx) %>%
#   separate(regionXdx2, into = c("BrainRegion","Dx"))
# 
# 
# 
# gene_sd_boxplot <- ggplot(gene_sd_long, aes(x = regionXdx, y = gene_sd, fill = BrainRegion)) +
#   geom_boxplot() +
#   theme_bw(base_size = 15) +
#   scale_fill_manual(values = region_colors)+
#   labs(x = "Region + Dx", y = "sd Gene expr_gene")
# 
# ggsave(here("case_control","plots","regionDx_var_boxplot.png"), width = 10)
# ggsave(here("case_control","plots","regionDx_var_boxplot.pdf"), width = 10)


#### Region Only ####
region <- levels(rse_gene$BrainRegion)
names(region) <- levels(rse_gene$BrainRegion)

region_sd <- map(expr, ~map_dfc(region, function(r){
  expr_gene_subset <- .x[,rse_gene$BrainRegion == r[[1]]]
  return(rowSds(expr_gene_subset))
}))

map(region_sd, head)

region_sd_long <- map2(region_sd, data, ~.x %>%
  add_column(transcript = rownames(.y)) %>%
  pivot_longer(!transcript, names_to = "region", values_to = "SD"))

walk2(region_sd_long, names(region_sd_long), function(sd, n){
  
  gene_sd_boxplot <- ggplot(sd, aes(x = region, y = SD, fill = region)) +
    geom_boxplot() +
    theme_bw(base_size = 15) +
    scale_fill_manual(values = region_colors)+
    labs(x = "Region", y = paste("Standard Deviation",n,"Expression"))+ 
    theme(legend.position = "none")
  
  fn = here("case_control","plots",paste0("region_var_boxplot-",tolower(n)))
  
  ggsave(gene_sd_boxplot, filename = paste0(fn, ".png"))
  ggsave(gene_sd_boxplot, filename = paste0(fn, ".pdf"))
  
})


# sgejobs::job_single('region_variability', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript region_variability.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-05-04 13:03:09 EDT"
# user  system elapsed 
# 37.290   2.305  44.019 
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
# date     2021-05-04                            
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# AnnotationDbi          1.52.0   2020-10-27 [2] Bioconductor                             
# askpass                1.1      2019-01-13 [2] CRAN (R 4.0.3)                           
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.1    2020-12-09 [1] CRAN (R 4.0.4)                           
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocFileCache          1.14.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.1   2021-04-16 [2] Bioconductor                             
# BiocParallel           1.24.1   2020-11-06 [2] Bioconductor                             
# biomaRt                2.46.3   2021-02-09 [2] Bioconductor                             
# Biostrings             2.58.0   2020-10-27 [2] Bioconductor                             
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.0.3)                           
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.0.3)                           
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.0.4)                           
# blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.3)                           
# broom                  0.7.6    2021-04-05 [2] CRAN (R 4.0.4)                           
# BSgenome               1.58.0   2020-10-27 [2] Bioconductor                             
# bumphunter             1.32.0   2020-10-27 [2] Bioconductor                             
# cachem                 1.0.4    2021-02-13 [2] CRAN (R 4.0.4)                           
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.3)                           
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
# dplyr                * 1.0.5    2021-03-05 [1] CRAN (R 4.0.4)                           
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.0.4)                           
# fansi                  0.4.2    2021-01-15 [2] CRAN (R 4.0.3)                           
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.0.4)                           
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.0.3)                           
# foreach                1.5.1    2020-10-15 [2] CRAN (R 4.0.3)                           
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.0.4)                           
# formatR                1.9      2021-04-14 [2] CRAN (R 4.0.4)                           
# Formula                1.2-4    2020-10-16 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.3)                           
# futile.logger        * 1.4.3    2016-07-10 [2] CRAN (R 4.0.3)                           
# futile.options         1.0.1    2018-04-20 [2] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.7   2021-04-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicAlignments      1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicFeatures        1.42.3   2021-04-01 [2] Bioconductor                             
# GenomicFiles           1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# GEOquery               2.58.0   2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.4.1    2021-04-23 [1] CRAN (R 4.0.4)                           
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
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.0.3)                           
# lambda.r               1.2.4    2019-09-18 [2] CRAN (R 4.0.3)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.4)                           
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.0.3)                           
# lifecycle              1.0.0    2021-02-15 [1] CRAN (R 4.0.4)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.3)                           
# lubridate              1.7.10   2021-02-26 [1] CRAN (R 4.0.4)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.4)                           
# MatrixGenerics       * 1.2.1    2021-01-30 [2] Bioconductor                             
# matrixStats          * 0.58.0   2021-01-29 [2] CRAN (R 4.0.3)                           
# memoise                2.0.0    2021-01-26 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# nnet                   7.3-15   2021-01-24 [3] CRAN (R 4.0.4)                           
# openssl                1.4.4    2021-04-30 [1] CRAN (R 4.0.4)                           
# pillar                 1.6.0    2021-04-13 [1] CRAN (R 4.0.4)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.3)                           
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.3)                           
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.3)                           
# ps                     1.6.0    2021-02-28 [1] CRAN (R 4.0.4)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# qvalue                 2.22.0   2020-10-27 [2] Bioconductor                             
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.0.3)                           
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.6    2021-01-15 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.0.4)                           
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.3)                           
# recount              * 1.16.1   2020-12-18 [2] Bioconductor                             
# rentrez                1.2.3    2020-11-10 [2] CRAN (R 4.0.3)                           
# reprex                 2.0.0    2021-04-02 [2] CRAN (R 4.0.4)                           
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                  0.4.11   2021-04-30 [1] CRAN (R 4.0.4)                           
# rngtools               1.5      2020-01-23 [2] CRAN (R 4.0.3)                           
# rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.0.4)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# Rsamtools              2.6.0    2020-10-27 [2] Bioconductor                             
# RSQLite                2.2.7    2021-04-22 [2] CRAN (R 4.0.4)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rtracklayer            1.50.0   2020-10-27 [2] Bioconductor                             
# rvest                  1.0.0    2021-03-09 [2] CRAN (R 4.0.4)                           
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.0.4)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# survival               3.2-9    2021-03-14 [3] CRAN (R 4.0.4)                           
# tibble               * 3.1.1    2021-04-18 [1] CRAN (R 4.0.4)                           
# tidyr                * 1.1.3    2021-03-03 [2] CRAN (R 4.0.4)                           
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.0.4)                           
# tidyverse            * 1.3.1    2021-04-15 [1] CRAN (R 4.0.4)                           
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.0.4)                           
# VariantAnnotation      1.36.0   2020-10-27 [2] Bioconductor                             
# vctrs                  0.3.8    2021-04-29 [1] CRAN (R 4.0.4)                           
# VennDiagram          * 1.6.20   2018-03-28 [2] CRAN (R 4.0.3)                           
# withr                  2.4.2    2021-04-18 [1] CRAN (R 4.0.4)                           
# xfun                   0.22     2021-03-11 [1] CRAN (R 4.0.4)                           
# XML                    3.99-0.6 2021-03-16 [2] CRAN (R 4.0.4)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
# **** Job ends ****
#   Tue May  4 13:03:10 EDT 2021
# 
