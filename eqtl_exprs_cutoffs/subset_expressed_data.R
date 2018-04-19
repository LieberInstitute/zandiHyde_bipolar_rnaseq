####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(segmented)
library(recount)

#################################  from get_expression_cutoffs.R
# Gene
# 2017-12-21 14:26:53 the suggested expression cutoff is 0.26 0.17
# Exon
# 2017-12-21 14:31:22 the suggested expression cutoff is 0.3 0.21
# Jxn
# 2017-12-21 14:35:11 the suggested expression cutoff is 0.21 0.33
# Tx
# 2017-12-21 14:37:48 the suggested expression cutoff is 0.4 0.24

### Final cutoffs used:
# Gene 0.25
# Exon 0.30
# Jxn 0.35
# Tx 0.40
############################################################


############################################################
##########  Amygdala

## load
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseGene_n511.rda")

regInd = which(colData(rse_gene)$BrainRegion == "Amygdala")

## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,regInd]

## exon
assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
exonIndex = rowMeans(assays(rse_exon)$rpkm) > 0.3
rse_exon = rse_exon[exonIndex,regInd]
mcols(rse_exon)$exonPos = paste0(seqnames(rse_exon),":",start(rse_exon),"-",end(rse_exon))

## junction
rowRanges(rse_jxn)$Length <- 100
assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
jxnIndex = rowMeans(assays(rse_jxn)$rp10m) > 0.35 & rowData(rse_jxn)$Class != "Novel"
rse_jxn = rse_jxn[jxnIndex,regInd]

## transcript
txIndex = rowMeans(assays(rse_tx)$tpm) > 0.4 
rse_tx = rse_tx[txIndex,regInd]

save(rse_gene,rse_exon,rse_jxn,rse_tx, file="eQTL_expressed_rse_amygdala.rda")

rse_gene_amyg = rse_gene
rse_exon_amyg = rse_exon
rse_jxn_amyg = rse_jxn
rse_tx_amyg = rse_tx
mcols(rse_jxn_amyg)$coord = paste0(ss(rownames(rse_jxn_amyg),"\\("), "(*)") ## to subset dlpfc later

############################################################
##########  sacc

## reload
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseGene_n511.rda")

regInd = which(colData(rse_gene)$BrainRegion == "sACC")

## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,regInd]

## exon
assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
exonIndex = rowMeans(assays(rse_exon)$rpkm) > 0.3
rse_exon = rse_exon[exonIndex,regInd]
mcols(rse_exon)$exonPos = paste0(seqnames(rse_exon),":",start(rse_exon),"-",end(rse_exon))

## junction
rowRanges(rse_jxn)$Length <- 100
assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
jxnIndex = rowMeans(assays(rse_jxn)$rp10m) > 0.35 & rowData(rse_jxn)$Class != "Novel"
rse_jxn = rse_jxn[jxnIndex,regInd]

## transcript
txIndex = rowMeans(assays(rse_tx)$tpm) > 0.4 
rse_tx = rse_tx[txIndex,regInd]

save(rse_gene,rse_exon,rse_jxn,rse_tx, file="eQTL_expressed_rse_sacc.rda")


############################################################
##########  DLPFC

## load
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")

## subset samples
subjInd = which(colData(rse_gene)$Race == "CAUC" & 
				colData(rse_gene)$Age > 13 &
				colData(rse_gene)$Dx %in% c("Bipolar","Control") )				

## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = which(rownames(rse_gene) %in% rownames(rse_gene_amyg) )
rse_gene = rse_gene[geneIndex,subjInd]
stopifnot(identical(rownames(rse_gene), rownames(rse_gene_amyg)))

## exon
assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
mcols(rse_exon)$exonPos = paste0(seqnames(rse_exon),":",start(rse_exon),"-",end(rse_exon))
exonIndex = which(rowData(rse_exon)$exonPos %in% rowData(rse_exon_amyg)$exonPos )
rse_exon = rse_exon[exonIndex,subjInd]
rse_exon = rse_exon[match(rowData(rse_exon_amyg)$exonPos, rowData(rse_exon)$exonPos),] ## match order of rows
stopifnot(identical(rowData(rse_exon)$exonPos, rowData(rse_exon_amyg)$exonPos))
mcols(rse_exon)$matchAmygStrand = ifelse(as.character(strand(rse_exon))==as.character(strand(rse_exon_amyg)),TRUE,FALSE)
table(mcols(rse_exon)$matchAmygStrand)   ## FALSE=3
mcols(rse_exon)$amygensemblID = mcols(rse_exon_amyg)$ensemblID
mcols(rse_exon)$amygSymbol = mcols(rse_exon_amyg)$Symbol

## junction
rowRanges(rse_jxn)$Length <- 100
assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
jxnIndex = which(rownames(rse_jxn) %in% mcols(rse_jxn_amyg)$coord )
rse_jxn = rse_jxn[jxnIndex,subjInd]

## transcript
txIndex = which(rownames(rse_tx) %in% rownames(rse_tx_amyg) ) 
rse_tx = rse_tx[txIndex,subjInd]
stopifnot(identical(rownames(rse_tx), rownames(rse_tx_amyg)))

save(rse_gene,rse_exon,rse_jxn,rse_tx, file="eQTL_expressed_rse_dlpfc.rda")











