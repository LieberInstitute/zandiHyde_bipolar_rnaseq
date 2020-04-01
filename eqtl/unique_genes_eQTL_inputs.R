
library(jaffelab)

#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

## Amyg and sACC
load("./raggr/mergedEqtl_output_sacc_raggr_4features.rda", verbose=TRUE)
sacc = allEqtl
jInd = which(!grepl("\\-", sacc$EnsemblGeneID) & sacc$Type=="Jxn")  ## cut ext from jxns
sacc$EnsemblGeneID[jInd] = ss(sacc$EnsemblGeneID[jInd], "\\.")

## cut after period for jxn Ensids

## DLPFC (has slightly different jxn numbers)
load("./raggr/mergedEqtl_output_dlpfc_raggr_4features.rda", verbose=TRUE)
d = allEqtl[which(allEqtl$Type == "Jxn"),]
d = d[!duplicated(d$EnsemblGeneID),]
d$EnsemblGeneID[!grepl("\\-", d$EnsemblGeneID)] = ss(d$EnsemblGeneID[!grepl("\\-", d$EnsemblGeneID)], "\\.")

## Create data frame
genes = data.frame(EnsemblGeneID = unique(c(sacc$EnsemblGeneID, d$EnsemblGeneID)), Symbol=NA)
genes$Symbol = sacc$Symbol[match(genes$EnsemblGeneID, sacc$EnsemblGeneID)]

genes$Gene881 = genes$EnsemblGeneID %in% sacc$EnsemblGeneID[sacc$Type=="Gene"]
genes$Exon881 = genes$EnsemblGeneID %in% sacc$EnsemblGeneID[sacc$Type=="Exon"]
genes$Jxn_sACCAmyg881 = genes$EnsemblGeneID %in% sacc$EnsemblGeneID[sacc$Type=="Jxn"]
genes$Jxn_DLPFC881 = genes$EnsemblGeneID %in% d$EnsemblGeneID
genes$Tx881 = genes$EnsemblGeneID %in% sacc$EnsemblGeneID[sacc$Type=="Tx"]

## How many unique genes were input to eQTL analyses
# table(genes$Gene881)
# FALSE  TRUE
 # 5252  4647

# table(genes$Exon881)
# FALSE  TRUE
 # 5454  4445

# table(genes$Jxn_sACCAmyg881)
# FALSE  TRUE
 # 6060  3839
 
# table(genes$Jxn_DLPFC881)
# FALSE  TRUE
 # 6099  3812

# table(genes$Tx881)
# FALSE  TRUE
 # 4880  5019



################
## Same thing but with inputs to smaller eQTL with 31 loci

sacc31 = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr_gwSignificant/all_raggr_31_snps_sacc_eqtls.csv",
	stringsAsFactors=FALSE)
jInd = which(!grepl("\\-", sacc31$EnsemblGeneID) & sacc31$Type=="Jxn")  ## cut ext from jxns
sacc31$EnsemblGeneID[jInd] = ss(sacc31$EnsemblGeneID[jInd], "\\.")

dlpfc31 = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr_gwSignificant/all_raggr_31_snps_dlpfc_eqtls.csv",
	stringsAsFactors=FALSE)
jInd = which(!grepl("\\-", dlpfc31$EnsemblGeneID))  ## cut ext from jxns
dlpfc31$EnsemblGeneID[jInd] = ss(dlpfc31$EnsemblGeneID[jInd], "\\.")

genes$Gene31 = genes$EnsemblGeneID %in% sacc31$EnsemblGeneID[sacc31$Type=="Gene"]
genes$Exon31 = genes$EnsemblGeneID %in% sacc31$EnsemblGeneID[sacc31$Type=="Exon"]
genes$Jxn_sACCAmyg31 = genes$EnsemblGeneID %in% sacc31$EnsemblGeneID[sacc31$Type=="Jxn"]
genes$Jxn_DLPFC31 = genes$EnsemblGeneID %in% dlpfc31$EnsemblGeneID
genes$Tx31 = genes$EnsemblGeneID %in% sacc31$EnsemblGeneID[sacc31$Type=="Tx"]


## How many unique genes were input to eQTL analyses
# table(genes$Gene31)
# FALSE  TRUE
 # 9454   457

# table(genes$Exon31)
# FALSE  TRUE
 # 9467   444

# table(genes$Jxn_sACCAmyg31)
# FALSE  TRUE
 # 9526   385

# table(genes$Jxn_DLPFC31)
# FALSE  TRUE
 # 5761   605

# table(genes$Tx31)
# FALSE  TRUE
 # 9406   505


 
 
# write.csv(genes, file="raggr_eQTL_inputs_unique_genes.csv")


################
## list index snp within 500 kb
library(SummarizedExperiment)

geneMap = genes[,1:2]
geneMap$id1 = ss(ss(as.character(geneMap$EnsemblGeneID),"-",1),"\\.")
geneMap$id2 = ss(ss(as.character(geneMap$EnsemblGeneID),"-",2),"\\.")
geneMap$id2[is.na(geneMap$id2)] = geneMap$id1[is.na(geneMap$id2)]


# gene coords
load("../preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata")
gMap = as.data.frame(rowRanges(rse_gene))
gMap = gMap[which(gMap$ensemblID %in% c(geneMap$id1,geneMap$id2)),1:10]

geneMap$chr = gMap$seqnames[match(geneMap$id1, gMap$ensemblID)]
geneMap$start1 = gMap$start[match(geneMap$id1, gMap$ensemblID)]
geneMap$start2 = gMap$start[match(geneMap$id2, gMap$ensemblID)]
geneMap$end1 = gMap$end[match(geneMap$id1, gMap$ensemblID)]
geneMap$end2 = gMap$end[match(geneMap$id2, gMap$ensemblID)]

geneMap$start = apply(geneMap[,6:7],1, FUN=min)
geneMap$end = apply(geneMap[,8:9],1, FUN=max)

geneMap = geneMap[,-c(6:9)]

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns=TRUE)



# risk snp coords
snpMap = read.csv("snpMap_10777.csv")
snpMap$IndexSNP = ss(as.character(snpMap$IndexSNP),":")
snpMapGR =  makeGRangesFromDataFrame(snpMap, seqnames.field="chr_hg38", 
				start.field="pos_hg38", end.field="pos_hg38",
				ignore.strand=TRUE, keep.extra.columns=TRUE)

				
				
## snps within 500 kb of each gene
f = findOverlaps(geneMapGR, snpMapGR,
	maxgap=5e5L, minoverlap=0L,
	type=c("any"),
	select=c("all"),
	ignore.strand=TRUE)
	
	
for (i in 1:nrow(geneMap)) {
	subjInd = subjectHits(f)[which(queryHits(f)==i)]
	snpMapSub = snpMap[subjInd,]
	inds = as.character(unique(snpMapSub$IndexSNP))
	if (length(inds)>0) {inds = paste(inds, collapse=",")
	} else { (inds=NA) }
	# coords = as.character(unique(snpMapSub$IndexSNP_hg38POS))
	# if (length(coords)>0) {coords = paste(coords, collapse=",")
	# } else { (coords=NA) }
	
	geneMap$indexSNP[i] = inds
	# geneMap$indexSNP_hg38POS[i] = coords
}




genes2 = genes
genes2$IndexSNP = geneMap$indexSNP


write.csv(genes2, file="raggr_eQTL_inputs_unique_genes_with_IndexSNP.csv")


















