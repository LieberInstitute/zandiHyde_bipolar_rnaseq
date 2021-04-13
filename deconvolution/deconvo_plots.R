
library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(reshape2)
library(patchwork)
library(tidyverse)
library(broom)
library(viridis)
library(sessioninfo)
library(patchwork)

## Load colors and plotting functions
cell_levels <- c("Astro","Micro","Oligo","OPC","Excit","Inhib")
cell_colors <- brewer.pal(n = length(cell_levels), name = "Set1")
names(cell_colors) <- cell_levels

## Venn Diagram colors
region_colors <- list(Amyg = "#ffff1f",
                      sACC = "#8eb0f6")
region_colors <- toupper(region_colors)

## Load MuSiC results & res data
load(here("deconvolution","est_prop_Bisque.Rdata"),verbose = TRUE)
load(here("data", "zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/sce_filtered.Rdata", verbose = TRUE)

long_prop <- est_prop_bisque$Est.prop.long %>%
    as_tibble() %>%
    mutate(cell_cat = case_when(grepl("Excit|Inhib",cell_type) ~ "Neuron",
                                TRUE ~ "Glia"))

long_prop$cell_type <- factor(long_prop$cell_type, levels = cell_levels)
levels(long_prop$cell_type)

pd <- as.data.frame(colData(rse_gene))
table(pd$PrimaryDx)
# Bipolar Control   Other 
# 245     264       2 
pd$PrimaryDx[pd$PrimaryDx == "Other"] <- "Bipolar"

pd2 <- pd[,c("Sex","PrimaryDx","BrainRegion")]%>%
  rownames_to_column("sample")

long_prop_pd <- left_join(long_prop, pd2, by = "sample")

## Boxplots
region_boxplot <- long_prop_pd %>%
  ggplot(aes(x = cell_type, y = prop, fill = BrainRegion)) +
  geom_boxplot() +
  scale_fill_manual(values = region_colors) +
  labs(x = "Cell Type", y = "Proportion", fill ='Brain Region') +
  theme_bw(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        base_size = 15)

ggsave(region_boxplot, filename = here("deconvolution","plots","cellType_boxplots_region.pdf"))

dx_boxPlot_all <- long_prop_pd %>%
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  facet_wrap(~BrainRegion)+
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    base_size = 15)

ggsave(dx_boxPlot_all, filename = here("deconvolution","plots","cellType_boxplots_dx.pdf"), width = 10)

#### composition barplot ####
sce_pd <- as.data.frame(colData(sce_pan))
mean_prop_sce <- sce_pd %>%
  select(cell_type = cellType.Broad, BrainRegion = region) %>%
  group_by(cell_type, BrainRegion) %>%
  summarise(n_cells = n()) %>%
  group_by(BrainRegion) %>%
  mutate(region_n_cells = sum(n_cells),
         mean_prop = n_cells/region_n_cells,
         data = "Single Cell Reference",
         BrainRegion = case_when(BrainRegion == "sacc" ~"sACC",
                                 BrainRegion == "amy"~ "Amygdala",
                                 TRUE ~ BrainRegion)) %>%
  filter(BrainRegion %in% c("sACC","Amygdala"))

mean_est_prop <- long_prop_pd %>%
  group_by(cell_type, BrainRegion) %>%
  summarize(mean_prop = mean(prop)) %>%
  mutate(data = "Bulk RNA-seq")
  
mean_prop <- rbind(mean_prop_sce, mean_est_prop) %>%
  arrange(cell_type) %>%
  group_by(data, BrainRegion) %>%
  mutate(anno_y = (1 - cumsum(mean_prop)) + mean_prop *.5)

comp_barplot <- mean_prop %>% ggplot(aes(x = BrainRegion, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~data) +
  scale_fill_manual(values = cell_colors)+
  geom_text(aes(y = anno_y, label = format(round(mean_prop,3),3)))+
  labs(x = "Brain Region", y = "Mean Proportion", fill ='Cell Type') +
  theme_bw(base_size = 15)

ggsave(comp_barplot, filename = here("deconvolution","plots","composition_barplot.pdf"))

#### Cor with qSV ####
load(here("case_control","qSV_mat.Rdata"), verbose = TRUE)
dim(qSV_mat)

qSV_long <- melt(qSV_mat) %>% rename(sample = Var1, qSV = Var2, qSV_value = value)

## Bind with qSV table
est_prop_qsv <- left_join(long_prop_pd, qSV_long, by = "sample")
levels(est_prop_qsv$cell_type)

#### Calculate p-values ####
prop_qSV_fit <- est_prop_qsv %>% group_by(cell_type, qSV, BrainRegion) %>%
  do(fitQSV = tidy(lm(prop ~ qSV_value-1, data = .))) %>% 
  unnest(fitQSV) %>%
  mutate(p.bonf = p.adjust(p.value, "bonf"),
         p.bonf.sig = p.bonf < 0.05,
         p.bonf.cat = cut(p.bonf, 
                          breaks = c(1,0.05, 0.01, 0.005, 0),
                          labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05")
                          ),
         p.fdr = p.adjust(p.value, "fdr"))

levels(prop_qSV_fit$p.bonf.cat)
prop_qSV_fit %>% count(BrainRegion, p.bonf.cat)
# BrainRegion p.bonf.cat     n
# <fct>       <fct>      <int>
# 1 Amygdala    <= 0.005      28
# 2 Amygdala    <= 0.01        3
# 3 Amygdala    <= 0.05        3
# 4 Amygdala    > 0.05        74
# 5 sACC        <= 0.005      31
# 6 sACC        <= 0.05        6
# 7 sACC        > 0.05        71

#### Save results ###
# save(prop_qSV_fit, file = here("deconvolution","prop_qSV_fit.Rdata"))

#### Tile plots ####
prop_qSV_fit1 <- prop_qSV_fit %>%
  mutate(log.p.bonf = -log10(p.bonf))

my_breaks <- c(0.05, 0.01, 0.005, 0)

sig_colors <- c(rev(viridis_pal(option = "magma")(3)),NA)
names(sig_colors) <- levels(prop_qSV_fit$p.bonf.cat)

tile_plot_val <-prop_qSV_fit1 %>%
  ggplot(aes(x = cell_type, y = qSV, fill = log.p.bonf)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = format(round(-log10(p.bonf),1), nsmall = 1),
                color = p.bonf.cat), size = 3, fontface = "bold", 
            show.legend = F)+
  scale_color_manual(values = sig_colors) +
  scale_fill_viridis(name = "-log10(p-value Bonf)", option = "magma", direction = -1) +
  facet_wrap(~BrainRegion)+ 
  labs(title ="p-values cell-type prop~qSV", x = 'Cell Type', color = "p-value Bonf\nsignificance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_bw(base_size = 15)

ggsave(plot = tile_plot_val,
       filename = here("deconvolution","plots","qSV_prop_fit_tileVal.pdf"), 
       width = 10)

#### Create scatter plots ####
# sig_colors2 <- c(brewer.pal(3, "Set1"),"black")
sig_colors2 <- c("#440154","#31688E", "#35B779","black")
names(sig_colors2) <- levels(prop_qSV_fit$p.bonf.cat)

regions <- list(sacc = "sACC", amyg = "Amygdala")

est_prop_qsv_fit <- left_join(est_prop_qsv, prop_qSV_fit)

scatter_plot <- map(regions, ~  est_prop_qsv_fit %>%
                      filter(BrainRegion == .x) %>%
                      ggplot(aes(x = qSV_value, y = prop, color = p.bonf.cat))+
                      geom_point(size = .4, alpha = .2) +
                      facet_grid(cell_type~qSV, scales = "free")+
                      theme_bw(base_size = 10)+
                      scale_color_manual(values = sig_colors2) +
                      theme(legend.text = element_text(size = 15)) +
                      guides(color = guide_legend(override.aes = list(size=5))) +
                      labs(title = .x)
)


ggsave(filename = here("deconvolution","plots", "qSV_cellType_scatter.png"),
       plot = scatter_plot$sacc/scatter_plot$amyg, width = 26, height = 12)

# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-01-11 11:12:17 EST"
# user  system elapsed 
# 42.305   2.392  46.937 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.0.3 Patched (2020-11-29 r79529)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-01-11                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.0    2020-11-02 [1] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0   2020-10-27 [2] Bioconductor                             
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.3)                           
# broom                * 0.7.3    2020-12-16 [2] CRAN (R 4.0.3)                           
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.3)                           
# cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)                           
# colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0    2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0   2020-10-27 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                * 1.0.2    2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.3)                           
# farver                 2.0.3    2020-01-16 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.0    2020-03-01 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2   2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1    2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0    2020-11-15 [1] CRAN (R 4.0.3)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0   2020-10-27 [1] Bioconductor                             
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jsonlite               1.7.1    2020-09-07 [1] CRAN (R 4.0.3)                           
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.0.3)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.3)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# lubridate              1.7.9.2  2020-11-13 [1] CRAN (R 4.0.3)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0    2020-10-27 [2] Bioconductor                             
# matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# patchwork            * 1.1.0    2020-11-09 [1] CRAN (R 4.0.3)                           
# pheatmap             * 1.0.12   2019-01-04 [2] CRAN (R 4.0.3)                           
# pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# ps                     1.4.0    2020-10-07 [1] CRAN (R 4.0.3)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)                           
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.3)                           
# reprex                 0.3.0    2019-05-16 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                * 0.4.9    2020-11-26 [1] CRAN (R 4.0.3)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rvest                  0.3.6    2020-07-25 [2] CRAN (R 4.0.3)                           
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-0    2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0   2020-10-27 [2] Bioconductor                             
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# tibble               * 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyr                * 1.1.2    2020-08-27 [2] CRAN (R 4.0.3)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.3)                           
# vctrs                  0.3.5    2020-11-17 [1] CRAN (R 4.0.3)                           
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
