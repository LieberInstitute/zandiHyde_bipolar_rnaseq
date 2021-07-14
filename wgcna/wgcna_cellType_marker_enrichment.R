library(SummarizedExperiment)
library(jaffelab)
library(here)
library(tidyverse)
library(WGCNA)
library(spatialLIBD)


load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
# filter ##
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)
marker_stats_filter <- marker_stats %>%
  filter(gene %in% rowData(rse_gene)$ensemblID) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  dplyr::slice(1:25)

marker_gene_list <- group_map(marker_stats_filter, ~pull(.x, gene))
names(marker_gene_list) <- levels(marker_stats_filter$cellType.target)

marker_stats_filter %>% summarise(max(rank_ratio))

load(here("wgcna","rdas","constructed_network_signed_bicor.rda"), verbose = TRUE)
# net_list
# net
# fNames

nrow(rse_gene) == length(net$colors)
# get colors
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
                      stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

colorDat
#      num          col Label numGenes
# 1      0         grey   ME0    13173
# 27     1    turquoise   ME1     2637
# 25     2         blue   ME2     2583
# 28     3        brown   ME3     1713
# 37     4       yellow   ME4     1440
# 30     5        green   ME5     1117
# 145    6          red   ME6      579
# 20     7        black   ME7      306
# 76     8         pink   ME8      224
# 40     9      magenta   ME9      220
# 39    10       purple  ME10      202
# 18    11  greenyellow  ME11      198
# 211   12          tan  ME12      139
# 285   13       salmon  ME13      102
# 1060  14         cyan  ME14       84
# 1172  15 midnightblue  ME15       82
# 188   16    lightcyan  ME16       78
# 443   17       grey60  ME17       69
# 26    18   lightgreen  ME18       66
# 153   19  lightyellow  ME19       65
# 654   20    royalblue  ME20       59

table(net$colorsLab == "grey",rowData(rse_gene)$ensemblID %in% marker_gene_list$Astro)
table(net$colorsLab == "red",rowData(rse_gene)$ensemblID %in% marker_gene_list$Micro)
#       FALSE  TRUE
# FALSE 24556     1
# TRUE    555    24
table(net$colorsLab == "pink",rowData(rse_gene)$ensemblID %in% marker_gene_list$Micro)
#       FALSE  TRUE
# FALSE 24887    25
# TRUE    224     0

table(net$colorsLab == "yellow",rowData(rse_gene)$ensemblID %in% marker_gene_list$Oligo)
#       FALSE  TRUE
# FALSE 23696     0
# TRUE   1415    25

chisq.test(table(net$colorsLab == "grey",rowData(rse_gene)$ensemblID %in% marker_gene_list$Oligo))
#       FALSE  TRUE
# FALSE 11938    25
# TRUE  13173     0
chisq.test(table(net$colorsLab == "grey",rowData(rse_gene)$ensemblID %in% marker_gene_list$Oligo))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(net$colorsLab == "grey", rowData(rse_gene)$ensemblID %in%     marker_gene_list$Oligo)
# X-squared = 25.493, df = 1, p-value = 4.44e-07
fisher.test(table(net$colorsLab == "grey",rowData(rse_gene)$ensemblID %in% marker_gene_list$Oligo))
# p-value = 8.568e-09

deModEnrich_markers <- map(marker_gene_list, function(markers) {
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net$colorsLab == cc,rowData(rse_gene)$ensemblID %in% markers)
    c(fisher.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  return(deModEnrich)
})

deModEnrich_table <- do.call("cbind", deModEnrich_markers)

deModEnrich_long <- deModEnrich_table %>% 
  rownames_to_column("ID") %>%
  pivot_longer(!ID, names_to = "stat") %>%
  separate(stat, into = c("test", "stat")) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  select(OR, Pval = Pvalue, ID, test) %>%
  filter(ID != "grey") %>%
  mutate(fdr = p.adjust(Pval, "fdr"),
         Pval.bonf = p.adjust(Pval, "bonferroni")) %>%
  as.data.frame()

deModEnrich_long %>% filter(fdr < 0.05) %>% arrange(fdr)
#            OR         Pval        ID  test          fdr    Pval.bonf
# 1 1061.881081 7.597905e-39       red Micro 1.367623e-36 1.367623e-36
# 2         Inf 7.342561e-32    yellow Oligo 6.608305e-30 1.321661e-29
# 3  108.177460 4.606883e-21     black Astro 2.764130e-19 8.292389e-19
# 4   22.602506 3.652416e-13      blue Inhib 1.643587e-11 6.574349e-11
# 5   66.380332 4.627276e-13   magenta Mural 1.665819e-11 8.329096e-11
# 6   18.241746 1.018357e-11 turquoise Excit 3.055071e-10 1.833043e-09
# 7  116.004630 3.070244e-09 royalblue  Endo 7.894913e-08 5.526439e-07
# 8   15.643486 1.319499e-03   magenta  Endo 2.968873e-02 2.375099e-01
# 9    8.127867 2.399039e-03       red Tcell 4.798079e-02 4.318271e-01

png(here("wgcna","gene_set_enrichment_cellType_markers.png"),
    width = 800, height = 600)
gene_set_enrichment_plot(deModEnrich_long)
title("OR: Top 25 Cell Type Markers & WGCNA Module")
dev.off()
