###
library(SummarizedExperiment)
library(clusterProfiler)

## load data
load("data/zandiHypde_bipolar_rseGene_n511.rda")

## read gene list
dat=  read.csv("cigenelist-yl20191217.csv",as.is=TRUE)

## get gene IDs
map = rowData(rse_gene)

eids = map$EntrezID[match(dat$ensemblid, map$ensemblID)]
eids = as.character(eids[!is.na(eids)])

univ = map$EntrezID
univ = as.character(univ[!is.na(univ)])

## run GO
go = enrichGO(eids, OrgDb = "org.Hs.eg.db", 
	ont="ALL",universe = univ, pvalueCutoff=1,qvalueCutoff=1)