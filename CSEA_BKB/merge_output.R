#

library(jaffelab)

## WGCNA
li_wgcna_files = list.files("csvs/Results", pattern = "_WGCNA_",full=TRUE)
li_wgcna_files = li_wgcna_files[grep("_Li", li_wgcna_files)]

names(li_wgcna_files) = ss(li_wgcna_files, "_", 3)

wgcnaList = lapply(li_wgcna_files, read.csv, as.is=TRUE,row.names=1,
	col.names = c("p05", "p01", "p001", "p0001"))

wgcnaTab = t(sapply(wgcnaList, function(x) x[,1]))
colnames(wgcnaTab) = rownames(wgcnaList[[1]])

options(width=100)
signif(wgcnaTab, 3)

## other analyses
li_files = list.files("csvs/Results", pattern = "_Li",full=TRUE)
li_files = li_files[!grepl("WGCNA", li_files)]

li_files_de = li_files[!grepl("_eQTL_", li_files)]
names(li_files_de) = ss(li_files_de, "/", 3)
names(li_files_de) = gsub("_Li.csv", "", names(li_files_de))

li_files_eqtl = li_files[grepl("_eQTL_", li_files)]
names(li_files_eqtl) = ss(li_files_eqtl, "/", 3)
names(li_files_eqtl) = gsub("_Li.csv", "", names(li_files_eqtl))

### DE
deList = lapply(li_files_de, read.csv, as.is=TRUE,row.names=1,
	col.names = c("p05", "p01", "p001", "p0001"))
deTab = t(sapply(deList, function(x) x[,1]))
colnames(deTab) = rownames(deList[[1]])

# eqtl
eqtlList = lapply(li_files_eqtl, read.csv, as.is=TRUE,row.names=1,
	col.names = c("p05", "p01", "p001", "p0001"))
eqtlTab = t(sapply(eqtlList, function(x) x[,1]))
colnames(eqtlTab) = rownames(eqtlList[[1]])