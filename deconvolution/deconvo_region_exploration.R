
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(xbioc)
library(BisqueRNA)
library(tidyverse)
library(reshape2)
library(DeconvoBuddies)
library(here)
library(sessioninfo)

#### Load Data ####
## Load rse_gene data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)
table(sce_pan$region)
# amy dlpfc   hpc   nac  sacc 
# 14039 11202 10139 19870 15323 

sce_pan <- sce_pan[, sce_pan$region %in% c("amy", "sacc")]
dim(sce_pan)
# [1] 23041 29362
length(unique(sce_pan$donor))
# [1] 5

table(sce_pan$region, sce_pan$cellType.Broad)
# Astro Endo Micro Mural Oligo  OPC Tcell Excit Inhib
# amy   1638   31  1201    39  6080 1459    31   443  3117
# sacc   907    0   784     0  4584  911     0  4163  3974
## missing Macro cell type

## Compute filtered region markers
markers_mean_ratio <- get_mean_ratio2(sce_pan, cellType_col = "cellType.Broad", assay_name = "logcounts", add_symbol = TRUE)

markers_stats_region <- markers_mean_ratio  %>%
  filter(gene %in% rownames(rse_gene)) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  slice(1:25)

marker_genes_region <- markers_stats_region$gene

## load other marker data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.Rdata", verbose = TRUE)
marker_stats_filter <- marker_stats %>%
  filter(gene %in% rownames(rse_gene)) %>%
  arrange(rank_ratio) %>%
  group_by(cellType.target) %>%
  slice(1:25)

length(intersect(markers_stats_region$gene, marker_stats_filter$gene))
# [1] 157

marker_stats_filter %>%
  mutate(region_marker = gene %in% markers_stats_region$gene) %>%
  group_by(cellType.target) %>%
  summarize(all_region = sum(!region_marker),
            common = sum(region_marker))

#### Run Bisque w/ this data ####
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene[marker_genes_region, ])$counts,
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce_pan[marker_genes_region,])$counts),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce_pan))[c("cellType.Broad", "cellType", "uniqueID","donor")]))

zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
exp_set_sce <- exp_set_sce[,zero_cell_filter]
# Exclude 5 cells

#### Run Bisque ####
est_prop_region <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "cellType.Broad",
                                        subject.names = "donor",
                                        use.overlap = FALSE)

est_prop_region$bulk.props <- t(est_prop_region$bulk.props)

est_prop_region$Est.prop.long <- melt(est_prop_region$bulk.props) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

#### Load full dataset bisque data ####
load(here("deconvolution","est_prop_Bisque.Rdata"), verbose = TRUE)

compare_props <- est_prop_region$Est.prop.long %>% 
  rename(region_prop = prop) %>%
  left_join(est_prop_bisque$Est.prop.long)

compare_props %>% 
  mutate(diff = abs(region_prop - prop)) %>%
  group_by(cell_type) %>% 
  summarise(mean_diff = mean(diff))

## Create color pallet
cell_colors <- create_cell_colors(pallet = "classic",
                                  cell_types = c("Astro", "Micro", "Oligo", "OPC", "Mural","Endo", "Tcell", "Inhib", "Excit"),
                                  preview = TRUE)

region_scatter <- ggplot(compare_props, aes(prop, region_prop, color = cell_type))+
  geom_point(size = 1) +
  geom_abline() +
  scale_color_manual(values = cell_colors) +
  labs(x = "Proportion w/ All Region Reference", y = "Proportion w/ sACC+Amyg Reference")

ggsave(region_scatter, filename = here("deconvolution", "plots", "region_specific_scatter.png"))

region_scatter_facet <- region_scatter +
  facet_wrap(~cell_type, scales = "free")+ 
  theme(aspect.ratio = 1)

ggsave(region_scatter_facet, filename = here("deconvolution", "plots", "region_specific_scatter_facet.png"), height = 10, width = 10)

compare_props %>% group_by(sample) %>%
  summarise(sum_region = sum(region_prop),
            sum = sum(prop)) %>%
  arrange(sum)
