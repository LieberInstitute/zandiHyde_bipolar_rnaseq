
library(SummarizedExperiment)
library(jaffelab)
library(matrixStats)
library(recount)
library(here)
library(sessioninfo)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)

source("leveneFast.R")
## load data
load(here("data","zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)

rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"),
                     levels = c("Control", "Bipolar"))

assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25
rse_gene = rse_gene[geneIndex,]

dim(rse_gene)
# [1] 25136   511

## add ancestry
load(here("genotype_data","zandiHyde_bipolar_MDS_n511.rda"), verbose = TRUE)
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## modJoint
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                        data=colData(rse_gene))

modJoint <- modJoint[,c(1:3,14,4:13)]
colnames(modJoint)

##Run cleaningY
expr <- cleaningY(log2(assays(rse_gene)$rpkm + 1), mod = modJoint, P = 4)

#### Levene Test ####

groups = factor(paste0(rse_gene$Dx, "_", rse_gene$BrainRegion))
table(groups)
# Bipolar_Amygdala     Bipolar_sACC Control_Amygdala     Control_sACC
#              121              126              122              142

levene_test_regionXdx <- leveneFast(dat = expr, groups = groups)
levene_test_regionXdx$adj.pval <- p.adjust(levene_test_regionXdx$pval, method = "BH")
table(levene_test_regionXdx$adj.pval < 0.05)
# FALSE  TRUE
# 20909  4227

levene_test_dx <- leveneFast(dat = expr, groups = rse_gene$Dx)
levene_test_dx$adj.pval <- p.adjust(levene_test_dx$pval, method = "BH")
table(levene_test_dx$adj.pval < 0.05)
# FALSE  TRUE
# 24406   730

levene_test_region <- leveneFast(dat = expr, groups = rse_gene$BrainRegion)
levene_test_region$adj.pval <- p.adjust(levene_test_region$pval, method = "BH")
table(levene_test_region$adj.pval < 0.05)
# FALSE  TRUE
# 20503  4633

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(
    rownames(levene_test_regionXdx)[levene_test_regionXdx$adj.pval < 0.05],
    rownames(levene_test_dx)[levene_test_dx$adj.pval < 0.05],
    rownames(levene_test_region)[levene_test_region$adj.pval < 0.05]),
  category.names = c("Region + Dx", 'Dx',"Region"),

  # Circles
  lwd = 2,
  fill = myCol,

  filename = here("case_control","plots","levene_test.png"),
  output = TRUE
)

# load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)
# statOut_gene <- statOut[statOut$Type == "Gene",]
# dim(statOut_gene)
#
#
# table(levene_test_regionXdx$adj.pval < 0.05, statOut_gene$adj.P.Val_sACC < 0.05)
# #       FALSE  TRUE
# # FALSE 20480   429
# # TRUE   3990   237
#
# table(levene_test_regionXdx$adj.pval < 0.05, statOut_gene$adj.P.Val_Amyg < 0.05)
# #       FALSE  TRUE
# # FALSE 20844    65
# # TRUE   4188    39
#
# venn.diagram(
#   x = list(
#     rownames(levene_test_regionXdx)[levene_test_regionXdx$adj.pval < 0.05],
#     rownames(statOut_gene)[statOut_gene$adj.P.Val_Amyg < 0.05],
#     rownames(statOut_gene)[statOut_gene$adj.P.Val_sACC < 0.05]),
#   category.names = c("Reg. + Dx LT", "DE Amyg","DE sACC"),
#
#   # Circles
#   lwd = 2,
#   fill = myCol,
#
#   filename = here("case_control","levene_test_DE.png"),
#   output = TRUE
# )


#### Box Plots ####
# it’s going to be 4 boxes (region * Dx), with the rowSds() output from the
# cleaningY(log2(RPKM + 1)) keeping Dx and Region effects and removing everything else (except the intercept)

regionXdx <- cross2(levels(rse_gene$BrainRegion),levels(rse_gene$Dx))
names(regionXdx) <- map(regionXdx, ~paste(.x[[1]],.x[[2]]))

table(rse_gene$BrainRegion, rse_gene$Dx)

gene_sd <- map_dfc(regionXdx, function(rXd){
  expr_subset <- expr[,rse_gene$BrainRegion == rXd[[1]] & rse_gene$Dx == rXd[[2]]]
  return(rowSds(expr_subset))
})

gene_sd_long <- gene_sd %>%
  add_column(gene = rownames(expr)) %>%
  pivot_longer(!gene, names_to = "regionXdx", values_to = "gene_sd") %>%
  mutate(regionXdx2 = regionXdx) %>%
  separate(regionXdx2, into = c("BrainRegion","Dx"))

region_colors <- list(Amyg = "#ffff1f",
                      sACC = "#8eb0f6")
region_colors <- toupper(region_colors)

gene_sd_boxplot <- ggplot(gene_sd_long, aes(x = regionXdx, y = gene_sd, fill = BrainRegion)) +
  geom_boxplot() +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = region_colors)+
  labs(x = "Region + Dx", y = "sd Gene Expr")

ggsave(here("case_control","plots","region_var_boxplot.png"), width = 10)
ggsave(here("case_control","plots","region_var_boxplot.pdf"), width = 10)

# sgejobs::job_single('region_variability', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript region_variability.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
