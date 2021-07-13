
library(tidyverse)
library(EnhancedVolcano)
library(here)
library(sessioninfo)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

statOut_long <- as.data.frame(statOut) %>% 
  select(starts_with("logFC"),starts_with("adj.P.Val"),Type, Symbol) %>%
  rownames_to_column("mol") %>%
  pivot_longer(!c(mol,Type, Symbol)) %>% 
  separate(name,sep = "_", into = c("var","Region")) %>% 
  pivot_wider(names_from = "var", values_from = "value")

statOut_long$Type <- factor(statOut_long$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))
head(statOut_long)

statOut_long %>% group_by(Region, Type) %>% 
  filter(abs(logFC) >1, adj.P.Val < 0.05) %>%
  summarize(n = n(),
            n_symb = paste(unique(Symbol), collapse = ", "))
# # A tibble: 6 x 4
# # Groups:   Region [2]
# Region Type         n n_symb              
# <chr>  <fct>    <int> <chr>               
#   1 Amyg   Gene         1 SLC5A5              
# 2 Amyg   Exon        14 IL1RL1, NR4A2, NR4A1
# 3 Amyg   Junction     1 NR4A1               
# 4 sACC   Gene         1 SERPINA3            
# 5 sACC   Exon         2 SLC11A1, HAMP       
# 6 sACC   Junction     2 HAMP, SLC11A1    

statOut_long %>%
  group_by(Region, Type) %>%
  summarise(max_lfc = max(logFC),
          min_lfc = min(logFC))

en_volcano <-  EnhancedVolcano(statOut_long,
                               lab = statOut_long$Symbol,
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                               pCutoff = 0.05) +
  facet_grid(Region~Type) +
  xlim(-1.7, 1.7)+
  ylim(0, 5)+
  theme_bw(base_size = 10)

ggsave(en_volcano, filename = here("case_control","plots","volcano_en.png"), width = 20, height = 10)
ggsave(en_volcano, filename = here("case_control","plots","volcano_en.pdf"), width = 20, height = 10)

#### Simple Volcano ####

simple_volcano <- ggplot(statOut_long, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.05) + 
  facet_grid(Region~Type)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "lightskyblue3")+
  theme_bw(base_size = 10) +
  labs(x = "log2 Fold Change", y = "-log10 FDR")

ggsave(simple_volcano, filename = here("case_control","plots","volcano_simple.png"), width = 10, height = 5)
ggsave(simple_volcano, filename = here("case_control","plots","volcano_simple.pdf"), width = 10, height = 5)

# sgejobs::job_single('volcano_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript volcano_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-04-26 10:21:33 EDT"
# user  system elapsed 
# 298.989   4.222 304.644 
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
# date     2021-04-26                            
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package         * version  date       lib source        
# ash               1.0-15   2015-09-01 [2] CRAN (R 4.0.3)
# assertthat        0.2.1    2019-03-21 [2] CRAN (R 4.0.3)
# backports         1.2.1    2020-12-09 [1] CRAN (R 4.0.4)
# beeswarm          0.3.1    2021-03-07 [1] CRAN (R 4.0.4)
# BiocGenerics    * 0.36.1   2021-04-16 [2] Bioconductor  
# broom             0.7.6    2021-04-05 [2] CRAN (R 4.0.4)
# cellranger        1.1.0    2016-07-27 [2] CRAN (R 4.0.3)
# cli               2.4.0    2021-04-05 [1] CRAN (R 4.0.4)
# colorspace        2.0-0    2020-11-11 [2] CRAN (R 4.0.3)
# crayon            1.4.1    2021-02-08 [2] CRAN (R 4.0.3)
# DBI               1.1.1    2021-01-15 [2] CRAN (R 4.0.3)
# dbplyr            2.1.1    2021-04-06 [2] CRAN (R 4.0.4)
# digest            0.6.27   2020-10-24 [1] CRAN (R 4.0.3)
# dplyr           * 1.0.5    2021-03-05 [1] CRAN (R 4.0.4)
# ellipsis          0.3.1    2020-05-15 [2] CRAN (R 4.0.3)
# EnhancedVolcano * 1.8.0    2020-10-27 [1] Bioconductor  
# extrafont         0.17     2014-12-08 [1] CRAN (R 4.0.4)
# extrafontdb       1.0      2012-06-11 [1] CRAN (R 4.0.4)
# fansi             0.4.2    2021-01-15 [2] CRAN (R 4.0.3)
# farver            2.1.0    2021-02-28 [2] CRAN (R 4.0.4)
# forcats         * 0.5.1    2021-01-27 [2] CRAN (R 4.0.3)
# fs                1.5.0    2020-07-31 [1] CRAN (R 4.0.3)
# generics          0.1.0    2020-10-31 [2] CRAN (R 4.0.3)
# ggalt             0.4.0    2017-02-15 [1] CRAN (R 4.0.4)
# ggbeeswarm        0.6.0    2017-08-07 [1] CRAN (R 4.0.3)
# ggplot2         * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)
# ggrastr           0.2.3    2021-03-01 [1] CRAN (R 4.0.4)
# ggrepel         * 0.9.1    2021-01-15 [2] CRAN (R 4.0.3)
# glue              1.4.2    2020-08-27 [1] CRAN (R 4.0.3)
# gtable            0.3.0    2019-03-25 [2] CRAN (R 4.0.3)
# haven             2.4.0    2021-04-14 [1] CRAN (R 4.0.4)
# here            * 1.0.1    2020-12-13 [1] CRAN (R 4.0.3)
# hms               1.0.0    2021-01-13 [2] CRAN (R 4.0.3)
# httr              1.4.2    2020-07-20 [1] CRAN (R 4.0.3)
# jsonlite          1.7.2    2020-12-09 [1] CRAN (R 4.0.3)
# KernSmooth        2.23-18  2020-10-29 [3] CRAN (R 4.0.4)
# labeling          0.4.2    2020-10-20 [2] CRAN (R 4.0.3)
# lifecycle         1.0.0    2021-02-15 [1] CRAN (R 4.0.4)
# lubridate         1.7.10   2021-02-26 [1] CRAN (R 4.0.4)
# magrittr          2.0.1    2020-11-17 [2] CRAN (R 4.0.3)
# maps              3.3.0    2018-04-03 [2] CRAN (R 4.0.3)
# MASS              7.3-53.1 2021-02-12 [3] CRAN (R 4.0.4)
# modelr            0.1.8    2020-05-19 [2] CRAN (R 4.0.3)
# munsell           0.5.0    2018-06-12 [2] CRAN (R 4.0.3)
# pillar            1.6.0    2021-04-13 [1] CRAN (R 4.0.4)
# pkgconfig         2.0.3    2019-09-22 [2] CRAN (R 4.0.3)
# proj4             1.0-10.1 2021-01-26 [1] CRAN (R 4.0.4)
# ps                1.6.0    2021-02-28 [1] CRAN (R 4.0.4)
# purrr           * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)
# R6                2.5.0    2020-10-28 [1] CRAN (R 4.0.3)
# RColorBrewer      1.1-2    2014-12-07 [2] CRAN (R 4.0.3)
# Rcpp              1.0.6    2021-01-15 [1] CRAN (R 4.0.3)
# readr           * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)
# readxl            1.3.1    2019-03-13 [2] CRAN (R 4.0.3)
# reprex            2.0.0    2021-04-02 [2] CRAN (R 4.0.4)
# rlang             0.4.10   2020-12-30 [1] CRAN (R 4.0.4)
# rprojroot         2.0.2    2020-11-15 [2] CRAN (R 4.0.3)
# rstudioapi        0.13     2020-11-12 [2] CRAN (R 4.0.3)
# Rttf2pt1          1.3.8    2020-01-10 [1] CRAN (R 4.0.4)
# rvest             1.0.0    2021-03-09 [2] CRAN (R 4.0.4)
# S4Vectors       * 0.28.1   2020-12-09 [2] Bioconductor  
# scales            1.1.1    2020-05-11 [2] CRAN (R 4.0.3)
# sessioninfo     * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)
# stringi           1.5.3    2020-09-09 [1] CRAN (R 4.0.3)
# stringr         * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)
# tibble          * 3.1.1    2021-04-18 [1] CRAN (R 4.0.4)
# tidyr           * 1.1.3    2021-03-03 [2] CRAN (R 4.0.4)
# tidyselect        1.1.0    2020-05-11 [2] CRAN (R 4.0.3)
# tidyverse       * 1.3.1    2021-04-15 [1] CRAN (R 4.0.4)
# utf8              1.2.1    2021-03-12 [2] CRAN (R 4.0.4)
# vctrs             0.3.7    2021-03-29 [1] CRAN (R 4.0.4)
# vipor             0.4.5    2017-03-22 [1] CRAN (R 4.0.3)
# withr             2.4.2    2021-04-18 [1] CRAN (R 4.0.4)
# xml2              1.3.2    2020-04-23 [2] CRAN (R 4.0.3)
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
# 
