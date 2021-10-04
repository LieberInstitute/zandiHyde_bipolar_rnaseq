library(SingleCellExperiment)
library(scater)
library(here)
library(ggplot2)
library(DeconvoBuddies)

## Use V1 sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.Rdata", verbose = TRUE)

## Create color pallet
cell_colors <- create_cell_colors(pallet = "classic",
                                  cell_types = c("Astro", "Micro", "Oligo", "OPC", "Mural","Endo", "Tcell", "Inhib", "Excit"),
                                  preview = TRUE)


sce_pan <- sce_pan[,sce_pan$region %in% c("amy","sacc")]

genes_of_intrest <- c("SCN2A", "GRIN2A")

genes_of_intrest %in%  rowData(sce_pan)$gene_name

rownames(sce_pan) <- rowData(sce_pan)$gene_name
sce_pan$ctXregion <- paste(sce_pan$cellType.Broad, sce_pan$region)
table(sce_pan$ctXregion)

pe <- plotExpression(sce_pan,
               exprs_values = "logcounts",
               features = genes_of_intrest,
               x = "ctXregion",
               colour_by = "cellType.Broad",
               point_alpha = 0.5,
               point_size = 0.2,
               ncol = 1,
               add_legend = F) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3) +
  scale_color_manual(values = cell_colors, guide="none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Cell Type + Brain Region")
  
ggsave(pe, filename = here("deconvolution","plots", paste0("sce_violin_", paste(genes_of_intrest, collapse = "-"), ".png")))
       
