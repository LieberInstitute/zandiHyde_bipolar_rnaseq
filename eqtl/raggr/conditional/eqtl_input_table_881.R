####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")

## dlp snps
# load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
# snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## load SNP data
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")
s = snpspos


##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gene))[,c(1:3,10)]
posGene$name = rownames(posGene)
posGene = posGene[,c(5,1:4)]

##### exon level 
posExon = as.data.frame(rowRanges(rse_exon))[,c(1:3,7,10)]
posExon$name = rownames(posExon)
posExon = posExon[,c(6,1:5)]

##### junction level 
posJxn = as.data.frame(rowRanges(rse_jxn))[,c(1:3,9,11,18,19)]
posJxn$name = rownames(posJxn)
posJxn = posJxn[,c(8,1:7)]

##### transcript level 
posTx = as.data.frame(rowRanges(rse_tx))[,c(1:3,10,13)]
posTx$name = rownames(posTx)
posTx = posTx[,c(6,1:5)]




s$downstream = s$pos-5e5
s$upstream = s$pos+5e5
s$downstream[s$downstream < 0] = 0

posGene$inRange = FALSE

for (i in 1:nrow(posGene)) {

	f_chr = as.character(posGene$seqnames[i])
	f_start = posGene$start[i]
	f_end = posGene$end[i]

	num = sum(f_chr == s$chr & 											## same chr
			( (f_start>s$downstream & f_start<s$upstream) | 	# starts inside window OR
			  (f_end>s$downstream & f_end<s$upstream) |			# ends inside window OR
			  (f_start<s$downstream & f_end>s$upstream) ) )		# covers whole window
			  
	if (num>0) { posGene$inRange[i] = TRUE }
}

posExon$inRange = FALSE

for (i in 1:nrow(posExon)) {

	f_chr = as.character(posExon$seqnames[i])
	f_start = posExon$start[i]
	f_end = posExon$end[i]

	num = sum(f_chr == s$chr & 											## same chr
			( (f_start>s$downstream & f_start<s$upstream) | 	# starts inside window OR
			  (f_end>s$downstream & f_end<s$upstream) |			# ends inside window OR
			  (f_start<s$downstream & f_end>s$upstream) ) )		# covers whole window
			  
	if (num>0) { posExon$inRange[i] = TRUE }
}



posJxn$inRange = FALSE

for (i in 1:nrow(posJxn)) {

	f_chr = as.character(posJxn$seqnames[i])
	f_start = posJxn$start[i]
	f_end = posJxn$end[i]

	num = sum(f_chr == s$chr & 											## same chr
			( (f_start>s$downstream & f_start<s$upstream) | 	# starts inside window OR
			  (f_end>s$downstream & f_end<s$upstream) |			# ends inside window OR
			  (f_start<s$downstream & f_end>s$upstream) ) )		# covers whole window
			  
	if (num>0) { posJxn$inRange[i] = TRUE }
}



posTx$inRange = FALSE

for (i in 1:nrow(posTx)) {

	f_chr = as.character(posTx$seqnames[i])
	f_start = posTx$start[i]
	f_end = posTx$end[i]

	num = sum(f_chr == s$chr & 											## same chr
			( (f_start>s$downstream & f_start<s$upstream) | 	# starts inside window OR
			  (f_end>s$downstream & f_end<s$upstream) |			# ends inside window OR
			  (f_start<s$downstream & f_end>s$upstream) ) )		# covers whole window
			  
	if (num>0) { posTx$inRange[i] = TRUE }
}





table(posGene$inRange)

table(posTx$inRange)
length(unique(posTx$gene_id[which(posTx$inRange)]))

table(posExon$inRange)
length(unique(posExon$gencodeID[which(posExon$inRange)]))

table(posJxn$inRange)
length(unique(posJxn$newGeneID[which(posJxn$inRange)]))


# > table(posGene$inRange)
# FALSE  TRUE
# 20489  4647
# >
# > table(posTx$inRange)
# FALSE  TRUE
# 58780 14434
# > length(unique(posTx$gene_id[which(posTx$inRange)]))
# [1] 5019
# >
# > table(posExon$inRange)
 # FALSE   TRUE
# 320229  76589
# > length(unique(posExon$gencodeID[which(posExon$inRange)]))
# [1] 4445
# >
# > table(posJxn$inRange)
 # FALSE   TRUE
# 217009  49188
# > length(unique(posJxn$newGeneID[which(posJxn$inRange)]))
# [1] 3839


### dlpfc
# > table(posJxn$inRange)
 # FALSE   TRUE
# 210751  48023
# > length(unique(posJxn$newGeneID[which(posJxn$inRange)]))
# [1] 3812














	  