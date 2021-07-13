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

# what if its not same markers as we used for deconvo?
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

deModEnrich_markers <- map(marker_gene_list, function(markers) {
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net$colorsLab == cc,rowData(rse_gene)$ensemblID %in% markers)
    c(chisq.test(tab)$p.value, getOR(tab))
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
  mutate(fdr = p.adjust(Pval, "fdr")) %>%
  as.data.frame()

deModEnrich_long %>% filter(fdr < 0.05) %>% arrange(fdr)
#              OR          Pval           ID  test           fdr
# 1  1.061881e+03 2.405140e-205          red Micro 4.545715e-203
# 2  1.081775e+02 4.229358e-128        black Astro 3.996743e-126
# 3           Inf  8.627184e-88       yellow Oligo  5.435126e-86
# 4  1.160046e+02  2.505465e-75    royalblue  Endo  1.183832e-73
# 5  6.638033e+01  8.395293e-71      magenta Mural  3.173421e-69
# 6  2.260251e+01  7.620750e-23         blue Inhib  2.400536e-21
# 7  1.824175e+01  1.284366e-19    turquoise Excit  3.467788e-18
# 8  0.000000e+00  4.440064e-07         grey Oligo  1.048965e-05
# 9  2.720761e+01  6.441992e-07 midnightblue   OPC  1.352818e-05
# 10 1.564349e+01  9.552242e-07      magenta  Endo  1.805374e-05
# 11 3.776635e-02  3.345374e-06         grey Micro  5.268964e-05
# 12 3.776635e-02  3.345374e-06         grey Excit  5.268964e-05
# 13 2.174870e+01  1.072296e-05       salmon Excit  1.558954e-04
# 14 8.127867e+00  9.602235e-05          red Tcell  1.296302e-03
# 15 1.727289e-01  5.681121e-04         grey Inhib  7.158213e-03
# 16 2.267429e-01  2.321124e-03         grey Astro  2.741828e-02

deModEnrich_long %>% filter(ID == "grey", test == "Oligo")
deModEnrich_long %>% filter(OR == 0 & fdr < .9)

png(here("wgcna","gene_set_enrichment_cellType_markers.png"),
    width = 800, height = 600)
gene_set_enrichment_plot(deModEnrich_long)
title("OR: Top 25 Cell Type Markers & WGCNA Module")
dev.off()
