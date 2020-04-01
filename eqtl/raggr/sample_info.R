### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
pda = colData(rse_gene)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
pds = colData(rse_gene)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pdd = colData(rse_gene)

pds$PrimaryDx[pds$PrimaryDx=="Other"] = "Bipolar"
pda$PrimaryDx[pda$PrimaryDx=="Other"] = "Bipolar"

table(pds$PrimaryDx)
table(pda$PrimaryDx)
table(pdd$Dx)

table(pds$PrimaryDx, pds$Sex)
table(pda$PrimaryDx, pda$Sex)
table(pdd$Dx, pdd$Sex)



summary(pds$AgeDeath[pds$PrimaryDx=="Bipolar"])
sd(pds$AgeDeath[pds$PrimaryDx=="Bipolar"])
summary(pds$AgeDeath[pds$PrimaryDx=="Control"])
sd(pds$AgeDeath[pds$PrimaryDx=="Control"])

summary(pda$AgeDeath[pda$PrimaryDx=="Bipolar"])
sd(pda$AgeDeath[pda$PrimaryDx=="Bipolar"])
summary(pda$AgeDeath[pda$PrimaryDx=="Control"])
sd(pda$AgeDeath[pda$PrimaryDx=="Control"])

summary(pdd$Age[pdd$Dx=="Bipolar"])
sd(pdd$Age[pdd$Dx=="Bipolar"])
summary(pdd$Age[pdd$Dx=="Control"])
sd(pdd$Age[pdd$Dx=="Control"])


## Overall
sex = c(as.character(pds$Sex),as.character(pda$Sex),as.character(pdd$Sex))
table(sex)
age = c(pds$AgeDeath,pda$AgeDeath,pdd$Age)
summary(age)
sd(age)


library(VennDiagram)
library(RColorBrewer)

v = venn.diagram(list(sACC = unique(pds$BrNum), Amygdala = unique(pda$BrNum), DLPFC = unique(pdd$BrNum)), 
	fill = c("skyblue3", "slateblue", "lightcyan2"), main="", main.pos = c(.5, .05), 
	cat.dist = c(.045,.06,.03), cat.cex = 1.9, cex=3, margin = .1, filename = NULL)	
pdf('venn_BrNums.pdf', useDingbats = FALSE)
palette(brewer.pal(8,"Dark2"))
    grid.draw(v)
dev.off()




## Separate cases and controls

## total unique
length(unique(c(pds$BrNum[pds$PrimaryDx=="Bipolar"],pda$BrNum[pda$PrimaryDx=="Bipolar"],
				pdd$BrNum[pdd$Dx=="Bipolar"])))				# 167
length(unique(c(pds$BrNum[pds$PrimaryDx=="Control"],pda$BrNum[pda$PrimaryDx=="Control"],
				pdd$BrNum[pdd$Dx=="Control"])))				# 213

### Overlaps
v1 = venn.diagram(list(sACC = unique(pds$BrNum[pds$PrimaryDx=="Bipolar"]), 
					Amygdala = unique(pda$BrNum[pda$PrimaryDx=="Bipolar"]), 
					DLPFC = unique(pdd$BrNum[pdd$Dx=="Bipolar"]) ), 
	fill = c("skyblue3", "slateblue", "lightcyan2"), main="", main.pos = c(.5, .05), 
	cat.dist = c(.045,.06,.03), cat.cex = 1.9, cex=3, margin = .1, filename = NULL)	
	
v2 = venn.diagram(list(sACC = unique(pds$BrNum[pds$PrimaryDx=="Control"]), 
					Amygdala = unique(pda$BrNum[pda$PrimaryDx=="Control"]), 
					DLPFC = unique(pdd$BrNum[pdd$Dx=="Control"]) ), 
	fill = c("skyblue3", "slateblue", "lightcyan2"), main="", main.pos = c(.5, .05), 
	cat.dist = c(.045,.06,.03), cat.cex = 1.9, cex=3, margin = .1, filename = NULL)	
	
pdf('venn_BrNums_cases_controls.pdf', useDingbats = FALSE)
palette(brewer.pal(8,"Dark2"))
    grid.draw(v1)
	grid.newpage()
	grid.draw(v2)
dev.off()







