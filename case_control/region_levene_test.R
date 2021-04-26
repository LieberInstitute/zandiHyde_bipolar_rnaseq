
library(SummarizedExperiment)
library(jaffelab)
library(matrixStats)
library(recount)
library(here)
library(sessioninfo)
library(VennDiagram)
library(RColorBrewer)

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

##fix coef
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
  
  filename = here("case_control","levene_test.png"),
  output = TRUE
)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)
statOut_gene <- statOut[statOut$Type == "Gene",]
dim(statOut_gene)


table(levene_test_regionXdx$adj.pval < 0.05, statOut_gene$adj.P.Val_sACC < 0.05)
#       FALSE  TRUE
# FALSE 20480   429
# TRUE   3990   237

table(levene_test_regionXdx$adj.pval < 0.05, statOut_gene$adj.P.Val_Amyg < 0.05)
#       FALSE  TRUE
# FALSE 20844    65
# TRUE   4188    39

venn.diagram(
  x = list(
    rownames(levene_test_regionXdx)[levene_test_regionXdx$adj.pval < 0.05],
    rownames(statOut_gene)[statOut_gene$adj.P.Val_Amyg < 0.05],
    rownames(statOut_gene)[statOut_gene$adj.P.Val_sACC < 0.05]),
  category.names = c("Reg. + Dx LT", "DE Amyg","DE sACC"),
  
  # Circles
  lwd = 2,
  fill = myCol,
  
  filename = here("case_control","levene_test_DE.png"),
  output = TRUE
)

# venn.diagram(
#   x = list(
#     rownames(levene_test_regionXdx)[levene_test_regionXdx$adj.pval < 0.05],
#     rownames(levene_test_region)[levene_test_region$adj.pval < 0.05],
#     rownames(statOut_gene)[statOut_gene$adj.P.Val_Amyg < 0.05],
#     rownames(statOut_gene)[statOut_gene$adj.P.Val_sACC < 0.05]),
#   category.names = c("Reg. + Dx LT","Region LT","DE Amyg","DE sACC"),
#   filename = here("case_control","levene_test_DE_region.png"),
#   output = TRUE
# )


# sgejobs::job_single('region_levene_test', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript region_levene_test.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
