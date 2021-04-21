
library(tidyverse)
library(spatialLIBD)
library(jaffelab)
library(here)


load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)
statOut <- statOut %>% 
  as.data.frame() %>% 
  filter(Type == "Gene") %>%
  mutate(ensemblID = ss(gencodeID,"\\."))

dim(statOut)

## marker data
load(here("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/marker_stats.Rdata"), verbose = TRUE)

marker_stats_50 <- marker_stats %>% 
  group_by(cellType.target) %>%
  filter(rank_ratio <= 50)

all(marker_stats_50$gene %in% statOut$ensemblID)
# [1] TRUE

tstats <- statOut[, grep("[f|t]_", colnames(statOut))]
fdrs <- statOut[, grep("adj.P.Val_", colnames(statOut))]
fdr_cut <- 0.05

geneList_present <- marker_stats_50 %>% 
  group_map(~pull(.x, gene))
names(geneList_present) <- levels(marker_stats_50$cellType.target)


source("gene_set_enrichment.R") 

gse <- gene_set_enrichment(gene_list = geneList_present, modeling_results = statOut)
gse %>% filter(Pval < 0.05)
#         OR         Pval   test Group direction fdr_cut
# 1 24.79029 3.073203e-14 t_sACC Micro      down    0.05

png("gene_set_erichment.png")
gene_set_enrichment_plot(gse)
dev.off()
