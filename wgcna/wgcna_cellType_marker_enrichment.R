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

n_genes <- list(`25` = 25, `50` = 50)

marker_stats_filter <- map(n_genes, ~marker_stats %>%
  filter(gene %in% rowData(rse_gene)$ensemblID) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  dplyr::slice(1:.x))

marker_gene_list <- map(marker_stats_filter,  function(ms){
  gene_list <- group_map(ms, ~pull(.x, gene))
  names(gene_list) <- levels(ms$cellType.target)
  return(gene_list)
})

map(marker_stats_filter, ~.x %>% summarise(max(rank_ratio)))

## Load wgcna data
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


table(net$colorsLab == "red",rowData(rse_gene)$ensemblID %in% marker_gene_list$`25`$Micro)
#       FALSE  TRUE
# FALSE 24556     1
# TRUE    555    24
table(net$colorsLab == "pink",rowData(rse_gene)$ensemblID %in% marker_gene_list$`25`$Micro)
#       FALSE  TRUE
# FALSE 24887    25
# TRUE    224     0

deModEnrich_markers <- map(marker_gene_list, ~map(.x, function(markers) {
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net$colorsLab == cc,rowData(rse_gene)$ensemblID %in% markers)
    c(fisher.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  return(deModEnrich)
}))

deModEnrich_long <- map(deModEnrich_markers, ~do.call("cbind", .x)  %>% 
  rownames_to_column("ID") %>%
  pivot_longer(!ID, names_to = "stat") %>%
  separate(stat, into = c("test", "stat")) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  select(OR, Pval = Pvalue, ID, test) %>%
  filter(ID != "grey") %>%
  mutate(fdr = p.adjust(Pval, "fdr"),
         Pval.bonf = p.adjust(Pval, "bonferroni")) %>%
  as.data.frame())

map(deModEnrich_long, ~.x %>% filter(fdr < 0.05) %>% arrange(fdr))
# $`25`
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
# 
# $`50`
#            OR         Pval         ID  test          fdr    Pval.bonf
# 1  281.356610 6.900566e-64        red Micro 1.242102e-61 1.242102e-61
# 2  834.690870 3.030843e-60     yellow Oligo 2.727759e-58 5.455517e-58
# 3   68.617203 2.308999e-29      black Astro 1.385399e-27 4.156198e-27
# 4   69.293317 6.321133e-25    magenta Mural 2.844510e-23 1.137804e-22
# 5   11.902242 2.555910e-16  turquoise Excit 9.201277e-15 4.600638e-14
# 6  109.914146 2.653745e-15  royalblue  Endo 7.961235e-14 4.776741e-13
# 7   10.347520 1.858113e-14       blue Inhib 4.778005e-13 3.344604e-12
# 8   10.771968 1.745819e-07        red Tcell 3.928092e-06 3.142474e-05
# 9   12.616438 8.215243e-05       pink Excit 1.478744e-03 1.478744e-02
# 10  12.853230 7.546134e-05    magenta  Endo 1.478744e-03 1.358304e-02
# 11  25.352584 3.102631e-04 lightgreen Excit 5.077033e-03 5.584736e-02
# 12  16.110251 1.107091e-03     salmon Excit 1.660636e-02 1.992763e-01
# 13   4.118167 1.511838e-03      green Inhib 2.093314e-02 2.721308e-01

## save all values
walk2(deModEnrich_long, names(deModEnrich_long), ~write.csv(.x, file = here("wgcna",paste0("wgcna_gene_set_enrichment_cellType_markers",.y,".csv")), row.names = FALSE))

walk2(deModEnrich_long, names(deModEnrich_long), function(g, n){
  name <- here("wgcna","plots",paste0("wgcna_gene_set_enrichment_cellType_markers",n))
  t <- paste("OR: Top",n,"Cell Type Markers & FDR < 0.05")
  
  png(paste0(name,".png"), width = 800)
  gene_set_enrichment_plot(g)
  title(t)
  dev.off()
  
  pdf(paste0(name, ".pdf"), width = 11)
  gene_set_enrichment_plot(g)
  title(t)
  dev.off()
  
})

# sgejobs::job_single('wgcna_cellType_marker_enrichment', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript wgcna_cellType_marker_enrichment.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

