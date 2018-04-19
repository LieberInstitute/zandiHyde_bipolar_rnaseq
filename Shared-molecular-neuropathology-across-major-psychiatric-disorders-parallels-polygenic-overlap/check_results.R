library(jaffelab)
library(limma)

## list tables
tabs = list.files("results/tables",full=TRUE, pattern=".csv")
names(tabs)= gsub(".csv", "", ss(tabs, "/",3))
names(tabs)= gsub("_092017", "", names(tabs))

## read in
tabList = lapply(tabs, read.csv,row.names=1,as.is=TRUE)

## get unique genes
g = unique(unlist(lapply(tabList, rownames)))
tabList = lapply(tabList, function(x) x[g,])
t(sapply(tabList,dim))

## split into microarray and RNAseq
tabListMicro = tabList[grep("^Micro", names(tabList))]
tabListSeq = tabList[grep("^RNA", names(tabList))]

## microarray only
fdrMatMicro = sapply(tabListMicro, "[[", "fdr")
fcMatMicro = sapply(tabListMicro, "[[", "beta")
colSums(fdrMatMicro < 0.05,na.rm=TRUE)

vennDiagram(vennCounts(fdrMatMicro[,-4] < 0.05))
vennDiagram(vennCounts(fdrMatMicro[,-4] < 0.01))

## rnaseq only
fdrMatSeq = lapply(tabListSeq, function(x) x[,grep("adj.P.Val", colnames(x))])
fdrMatSeq = do.call("cbind", fdrMatSeq)
pvalMatSeq = lapply(tabListSeq, function(x) x[,grep("P.Value", colnames(x))])
pvalMatSeq = do.call("cbind", pvalMatSeq)
tstatMatSeq = lapply(tabListSeq, function(x) x[,grep("t", colnames(x))])
tstatMatSeq = do.call("cbind", tstatMatSeq)

## via qSVA
fdrMatSeq_qSVA = fdrMatSeq[,grep("qsva", colnames(fdrMatSeq),ignore=TRUE)]
colnames(fdrMatSeq_qSVA) = gsub("RNAseq_", "", colnames(fdrMatSeq_qSVA))
colnames(fdrMatSeq_qSVA) = gsub("adj.P.Val", "", colnames(fdrMatSeq_qSVA))

vennDiagram(vennCounts(fdrMatSeq_qSVA < 0.05))
vennDiagram(vennCounts(fdrMatSeq_qSVA[,c(2,4,3,5)] < 0.05))
vennDiagram(vennCounts(fdrMatSeq_qSVA[,c(2,4,3,5)] < 0.05))

###################
## other modeling
###############
## just ASD
fdrMatSeq_ASD = fdrMatSeq[,1:5]
colnames(fdrMatSeq_ASD) = ss(colnames(fdrMatSeq_ASD), "\\.",2)
vennDiagram(vennCounts(fdrMatSeq_ASD < 0.05))

fdrMatSeq_SZ = fdrMatSeq[,grep(".SCZ.", colnames(fdrMatSeq))[seq(1,8,2)]]
colnames(fdrMatSeq_SZ) = gsub("RNAseq_", "", colnames(fdrMatSeq_SZ))
colnames(fdrMatSeq_SZ) = gsub("adj.P.Val", "", colnames(fdrMatSeq_SZ))
vennDiagram(vennCounts(fdrMatSeq_SZ < 0.05))

fdrMatSeq_BPD = fdrMatSeq[,grep(".SCZ.", colnames(fdrMatSeq))[seq(2,8,2)]]
colnames(fdrMatSeq_BPD) = gsub("RNAseq_", "", colnames(fdrMatSeq_BPD))
colnames(fdrMatSeq_BPD) = gsub("adj.P.Val", "", colnames(fdrMatSeq_BPD))
vennDiagram(vennCounts(fdrMatSeq_BPD < 0.05))

fdrMatSeq_Obs = fdrMatSeq[,c(1,7,8,11,12)]
colnames(fdrMatSeq_Obs) = gsub("RNAseq_", "", colnames(fdrMatSeq_Obs))
colnames(fdrMatSeq_Obs) = gsub("adj.P.Val", "", colnames(fdrMatSeq_Obs))
vennDiagram(vennCounts(fdrMatSeq_Obs < 0.05))
vennDiagram(vennCounts(fdrMatSeq_Obs[,-1] < 0.05))
vennDiagram(vennCounts(fdrMatSeq_Obs[,c(2,4,3,5)] < 0.05))

colnames(fdrMatSeq_qSVA) = gsub("RNAseq_", "", colnames(fdrMatSeq_qSVA))

colnames(fdrMatSeq_SZ) = ss(colnames(fdrMatSeq_SZ), "\\.",2)
vennDiagram(vennCounts(fdrMatSeq_ASD < 0.05))

# number of FDR
lapply(tabList,head,2)
