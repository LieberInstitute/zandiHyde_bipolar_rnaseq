library(readxl)
library(SummarizedExperiment)

## read in magma
magma = read_excel("../PGC_BP_Magma_table.xlsx", skip = 2)
magma = as.data.frame(magma)

## load DE
load("bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda")
geneStat = statOut[statOut$Type == "Gene",]

## line up
mm = match(magma$SYMBOL, geneStat$Symbol)
table(is.na(mm))
magma$SYMBOL[which(is.na(mm))[1:10]]

plot(geneStat$t_sACC[mm], magma$ZSTAT)
plot(-log10(geneStat$P.Value_sACC[mm]), -log10(magma$P_JOINT),
	xlab = "BPD DE- sACC", ylab = "Magma", main = "-log10 p-values")
plot(-log10(geneStat$P.Value_Amyg[mm]), -log10(magma$P_JOINT),
	xlab = "BPD DE- Amygdala", ylab = "Magma", main = "-log10 p-values")
