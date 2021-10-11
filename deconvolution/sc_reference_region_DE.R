library(SingleCellExperiment)
library(jaffelab)
library(edgeR)
library(scran)
library(scuttle)
library(here)
library(sessioninfo)
library(tidyverse)

## Use V1 sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)

## get marker genes
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
marker_stats_filter <- marker_stats %>%
  filter(gene %in% rownames(rse_gene)) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  slice(1:25) %>%
  select(gene, cellType.target, rank_ratio, Symbol)

table(marker_stats_filter$cellType.target)

## filter to regions of intrest
table(sce_pan$region)
table(sce_pan$region, sce_pan$donor)
sce_pan <- sce_pan[,sce_pan$region %in% c("amy","sacc")]

table(sce_pan$sampleID, sce_pan$cellType.Broad)

summed <- aggregateAcrossCells(sce_pan, 
                               id=colData(sce_pan)[,c("cellType.Broad", "sampleID")])

colData(summed)

#### preform DE w/ pseudobulk data ####
table(summed$cellType.Broad, summed$region)

summed.filt <- summed[,summed$ncells >= 10]
table(summed.filt$cellType.Broad, summed.filt$region)

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$cellType.Broad,
                            design= ~region + donor,
                            coef="regionsacc",
                            condition=summed.filt$region
)

names(de.results)
sum(de.results$Astro$FDR < 0.05, na.rm = TRUE)

length(de.results)

marker_stats_excit <- marker_stats_filter %>% filter(cellType.target == "Excit")

de.results.all <- as.data.frame(do.call("rbind", de.results)) %>%
  add_column(cellType.target = rep(names(de.results), each = nrow(de.results$Astro))) %>%
  rownames_to_column("gene") %>%
  mutate(gene = ss(gene,"\\.")) %>%
  left_join(marker_stats_filter)

de.results.all %>% 
  group_by(cellType.target) %>% 
  filter(FDR < 0.05) %>%
  count()

# cellType.target     n
# <chr>           <int>
#   1 Astro             351
# 2 Excit             879
# 3 Inhib             140
# 4 Micro              29
# 5 Oligo             531
# 6 OPC               183


write_csv(de.results.all, file = here("deconvolution","sc_reference_region_DE.csv"))

de_markers <- de.results.all %>%
  filter(!is.na(rank_ratio),FDR < 0.05) 
#               gene      logFC   logCPM         F       PValue          FDR cellType.target rank_ratio      Symbol
# 1  ENSG00000164089 -0.8711247 8.924453  42.82804 1.045172e-05 0.0017745201           Astro          1      ETNPPL
# 2  ENSG00000271904 -0.5097911 9.081093  15.37328 1.420176e-03 0.0479462964           Astro         15  AC091826.2
# 3  ENSG00000170324 -1.7568678 7.173895  60.44846 4.250527e-06 0.0009221282           Astro         20      FRMPD2
# 4  ENSG00000171885 -0.7105329 9.208903  16.55637 1.148303e-03 0.0423030318           Astro          9        AQP4
# 5  ENSG00000188039  0.6081967 7.431762  22.11105 3.018313e-04 0.0179191463           Astro         22        NWD1
# 6  ENSG00000239498  2.4723069 7.744938 124.67011 1.776554e-06 0.0004815502           Excit          5  AC019211.1
# 7  ENSG00000236451  2.3125868 7.294832 246.24022 1.341107e-08 0.0000151587           Excit          2  AC067956.1
# 8  ENSG00000236922  2.9876100 7.742432  66.34982 2.537313e-04 0.0113301998           Excit         21   LINC01378
# 9  ENSG00000146147  1.1729183 8.218120  46.60509 3.694294e-05 0.0031514768           Excit          8        MLIP
# 10 ENSG00000263878  2.1129302 6.034381 120.55469 4.578488e-07 0.0002070049           Excit         18  DLGAP1-AS4
# 11 ENSG00000141668  2.1988192 8.078708 116.19709 1.328376e-06 0.0004142014           Excit          3       CBLN2
# 12 ENSG00000235448  0.6747846 8.376845  18.92126 4.667003e-05 0.0029141960           Oligo         20 LURAP1L-AS1


write_csv(de_markers, file = here("deconvolution","sc_reference_region_DE_significant_markers.csv"))

de_markers %>% count(cellType.target)  
#   cellType.target n
# 1           Astro 5
# 2           Excit 6
# 3           Oligo 1

