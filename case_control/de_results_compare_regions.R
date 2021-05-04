
library(tidyverse)
library(here)
library(sessioninfo)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

head(statOut)
table(statOut$Type)

statOut <- as.data.frame(statOut) %>% 
  mutate(`FDR < 0.05` = case_when(adj.P.Val_sACC < 0.05 & adj.P.Val_Amyg < 0.05 ~"Both",
                                  adj.P.Val_sACC < 0.05 ~ "sACC",
                                  adj.P.Val_Amyg < 0.05 ~ "Amygdala",
                                  TRUE ~ "None"))

statOut$Type <- factor(statOut$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))
statOut$`FDR < 0.05` <- factor(statOut$`FDR < 0.05`, levels = c("Both", "Amygdala", "sACC", "None"))

table(statOut$`FDR < 0.05`)
# Both Amygdala     sACC     None 
# 67      253     1768   759277 

## spearman cor
region_cor <- statOut %>% 
  group_by(Type) %>%
  summarize(cor = cor(t_Amyg, t_sACC, method = "spearman")) %>%
  mutate(cor_anno = paste0("rho==", format(round(cor, 2), nsmall = 2)))

fdr_colors <- c(None = "gray", sACC = "skyblue3", Amygdala = "orange", Both = "purple")

t_stat_scatter <- ggplot(statOut, aes(x = t_Amyg, y = t_sACC)) +
  geom_point(aes(color = `FDR < 0.05`), size = 0.5, alpha = 0.5) +
  facet_wrap(~Type) +
  labs(x = "t-stat Amygdala", y = "t-stat sACC") +
  scale_color_manual(values = fdr_colors) +
  geom_text(data = region_cor, aes(x = -4.25, y = 6, label = cor_anno), parse = TRUE) +
  theme_bw(base_size = 15) +
  NULL

ggsave(t_stat_scatter, filename = here("case_control","plots","deStats_region_t_plots.png"))
ggsave(t_stat_scatter, filename = here("case_control","plots","deStats_region_t_plots.pdf"))

# sgejobs::job_single('de_results_compare_regions', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript de_results_compare_regions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-05-04 13:15:21 EDT"
# user  system elapsed 
# 89.654   1.539  91.743 
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
# package      * version date       lib source        
# assertthat     0.2.1   2019-03-21 [2] CRAN (R 4.0.3)
# backports      1.2.1   2020-12-09 [1] CRAN (R 4.0.4)
# BiocGenerics * 0.36.1  2021-04-16 [2] Bioconductor  
# broom          0.7.6   2021-04-05 [2] CRAN (R 4.0.4)
# cellranger     1.1.0   2016-07-27 [2] CRAN (R 4.0.3)
# cli            2.5.0   2021-04-26 [1] CRAN (R 4.0.4)
# colorspace     2.0-0   2020-11-11 [2] CRAN (R 4.0.3)
# crayon         1.4.1   2021-02-08 [2] CRAN (R 4.0.3)
# DBI            1.1.1   2021-01-15 [2] CRAN (R 4.0.3)
# dbplyr         2.1.1   2021-04-06 [2] CRAN (R 4.0.4)
# digest         0.6.27  2020-10-24 [1] CRAN (R 4.0.3)
# dplyr        * 1.0.5   2021-03-05 [1] CRAN (R 4.0.4)
# ellipsis       0.3.2   2021-04-29 [2] CRAN (R 4.0.4)
# fansi          0.4.2   2021-01-15 [2] CRAN (R 4.0.3)
# farver         2.1.0   2021-02-28 [2] CRAN (R 4.0.4)
# forcats      * 0.5.1   2021-01-27 [2] CRAN (R 4.0.3)
# fs             1.5.0   2020-07-31 [1] CRAN (R 4.0.3)
# generics       0.1.0   2020-10-31 [2] CRAN (R 4.0.3)
# ggplot2      * 3.3.3   2020-12-30 [2] CRAN (R 4.0.3)
# glue           1.4.2   2020-08-27 [1] CRAN (R 4.0.3)
# gtable         0.3.0   2019-03-25 [2] CRAN (R 4.0.3)
# haven          2.4.1   2021-04-23 [1] CRAN (R 4.0.4)
# here         * 1.0.1   2020-12-13 [1] CRAN (R 4.0.3)
# hms            1.0.0   2021-01-13 [2] CRAN (R 4.0.3)
# httr           1.4.2   2020-07-20 [1] CRAN (R 4.0.3)
# jsonlite       1.7.2   2020-12-09 [1] CRAN (R 4.0.3)
# labeling       0.4.2   2020-10-20 [2] CRAN (R 4.0.3)
# lifecycle      1.0.0   2021-02-15 [1] CRAN (R 4.0.4)
# lubridate      1.7.10  2021-02-26 [1] CRAN (R 4.0.4)
# magrittr       2.0.1   2020-11-17 [2] CRAN (R 4.0.3)
# modelr         0.1.8   2020-05-19 [2] CRAN (R 4.0.3)
# munsell        0.5.0   2018-06-12 [2] CRAN (R 4.0.3)
# pillar         1.6.0   2021-04-13 [1] CRAN (R 4.0.4)
# pkgconfig      2.0.3   2019-09-22 [2] CRAN (R 4.0.3)
# ps             1.6.0   2021-02-28 [1] CRAN (R 4.0.4)
# purrr        * 0.3.4   2020-04-17 [1] CRAN (R 4.0.3)
# R6             2.5.0   2020-10-28 [1] CRAN (R 4.0.3)
# Rcpp           1.0.6   2021-01-15 [1] CRAN (R 4.0.3)
# readr        * 1.4.0   2020-10-05 [2] CRAN (R 4.0.3)
# readxl         1.3.1   2019-03-13 [2] CRAN (R 4.0.3)
# reprex         2.0.0   2021-04-02 [2] CRAN (R 4.0.4)
# rlang          0.4.11  2021-04-30 [1] CRAN (R 4.0.4)
# rprojroot      2.0.2   2020-11-15 [2] CRAN (R 4.0.3)
# rstudioapi     0.13    2020-11-12 [2] CRAN (R 4.0.3)
# rvest          1.0.0   2021-03-09 [2] CRAN (R 4.0.4)
# S4Vectors    * 0.28.1  2020-12-09 [2] Bioconductor  
# scales         1.1.1   2020-05-11 [2] CRAN (R 4.0.3)
# sessioninfo  * 1.1.1   2018-11-05 [2] CRAN (R 4.0.3)
# stringi        1.5.3   2020-09-09 [1] CRAN (R 4.0.3)
# stringr      * 1.4.0   2019-02-10 [2] CRAN (R 4.0.3)
# tibble       * 3.1.1   2021-04-18 [1] CRAN (R 4.0.4)
# tidyr        * 1.1.3   2021-03-03 [2] CRAN (R 4.0.4)
# tidyselect     1.1.1   2021-04-30 [2] CRAN (R 4.0.4)
# tidyverse    * 1.3.1   2021-04-15 [1] CRAN (R 4.0.4)
# utf8           1.2.1   2021-03-12 [2] CRAN (R 4.0.4)
# vctrs          0.3.8   2021-04-29 [1] CRAN (R 4.0.4)
# withr          2.4.2   2021-04-18 [1] CRAN (R 4.0.4)
# xml2           1.3.2   2020-04-23 [2] CRAN (R 4.0.3)
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
