##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)

################
## load SNP data
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## risk loci from PGC paper
indexLoci = read.csv("../../PGC_risk_loci.csv", stringsAsFactors=FALSE) ## 881
indexIndex = which(snpMap$pos_hg19 %in% indexLoci$hg19POS)	# keep 456

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

riskLoci$genomewideSig1 = NA
for (i in 1:nrow(riskLoci)) {
	indexInd = which(indexLoci$hg19POS == riskLoci$hg19POS1[i])
	riskLoci$genomewideSig1[i] = indexLoci$genomewide[indexInd]
}
riskLoci = riskLoci[which(riskLoci$genomewideSig1==TRUE),]
riskLoci_full = riskLoci

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 1226 snps from snpMap
snpMap = snpMap[keepIndex,]
keepIndex = which(riskLoci$hg19POS2 %in% snpMap$pos_hg19)	# keep 1223 snps from riskLoci
riskLoci = riskLoci[keepIndex,]

snpMap$Status = ifelse(snpMap$pos_hg19 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status1 = ifelse(riskLoci$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status2 = ifelse(riskLoci$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")

## Also keep track of full list before dropping
riskLoci_full$hg19POS1 = paste0(riskLoci_full$SNP1_Chr, ":", riskLoci_full$SNP1_Pos) 
riskLoci_full$hg19POS2 = paste0(riskLoci_full$SNP2_Chr, ":", riskLoci_full$SNP2_Pos) 
riskLoci_full$SNP2_missing = "missing"
riskLoci_full$SNP2_missing[keepIndex] = "analyzed"
riskLoci_full$Status1 = ifelse(riskLoci_full$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci_full$Status2 = ifelse(riskLoci_full$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")


### numbers:
# 27 indexSNPs											# length(unique(riskLoci$SNP1_Name))
# 456 indexSNPs in our SNPs used in eQTLS				#  length(unique(riskLoci$hg19POS1)) 490?
# 13,592 rAggr riskLoci (including index SNPs)  		#  nrow(riskLoci_full)
# 1226 rAggr riskLoci in our SNPs used in eQTLS 		#  nrow(snpMap)


################
## load table
amyg = read.csv("raggr_genomewide31_snps_amyg_eqtls_fdr01.csv", row.names=1)
sacc = read.csv("raggr_genomewide31_snps_sacc_eqtls_fdr01.csv", row.names=1)
dlp = read.csv("raggr_genomewide31_snps_dlpfc_eqtls_fdr01.csv", row.names=1)

## unique SNPs
length(unique(sacc$hg19POS))
length(unique(amyg$hg19POS))
length(unique(dlp$hg19POS))

venn.diagram(list(Amygdala = amyg$hg19POS, sACC = sacc$hg19POS, DLPFC = dlp$hg19POS), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_SNP.png")

## unique index SNPs
length(unique(sacc$hg19POS[sacc$hg19POS %in% riskLoci$hg19POS1]))
length(unique(amyg$hg19POS[amyg$hg19POS %in% riskLoci$hg19POS1]))
length(unique(dlp$hg19POS[dlp$hg19POS %in% riskLoci$hg19POS1]))

venn.diagram(list(Amygdala = unique(amyg$hg19POS[amyg$hg19POS %in% riskLoci$hg19POS1]), 
				sACC = unique(sacc$hg19POS[sacc$hg19POS %in% riskLoci$hg19POS1]), 
				DLPFC = unique(dlp$hg19POS[dlp$hg19POS %in% riskLoci$hg19POS1])), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_index.png")

	
## unique features
tapply(sacc$gene, sacc$Type, function(x) length(unique(x)))
tapply(amyg$gene, amyg$Type, function(x) length(unique(x)))
tapply(dlp$gene, dlp$Type, function(x) length(unique(x)))

venn.diagram(list(Amygdala = amyg$gene, sACC = sacc$gene, DLPFC = dlp$gene), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_features.png")

## SNP-feature pairs
nrow(sacc)  ## 38609
nrow(amyg)   ## 22572
nrow(dlp)   ## 14669
table(sacc$Type)
table(amyg$Type)
table(dlp$Type)

venn.diagram(list(Amygdala = paste0(amyg$SNP, amyg$gene), 
				sACC = paste0(sacc$SNP, sacc$gene), 
				DLPFC = paste0(dlp$SNP, dlp$gene) ), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_snpfeatpairs.png")



## Unique symbols in SNP-feature pairs
length(unique(sacc$Symbol))
length(unique(amyg$Symbol))
length(unique(dlp$Symbol))
tapply(sacc$Symbol, sacc$Type, function(x) length(unique(x)))
tapply(amyg$Symbol, amyg$Type, function(x) length(unique(x)))
tapply(dlp$Symbol, dlp$Type, function(x) length(unique(x)))

venn.diagram(list(Amygdala = amyg$Symbol, sACC = sacc$Symbol, DLPFC = dlp$Symbol), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_symbol.png")

	
	
	
	
################################################################
###### Index SNP info ##########
################################################################

region = amyg

## note which proxy snps have a significant result
riskLoci$proxy_FDRsig = "na"
for (i in 1:nrow(riskLoci)) {
	pos = riskLoci$hg19POS2[i]
	sig = which(region$hg19POS == pos)
	riskLoci$proxy_FDRsig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
}	

## Is SNP index snp or proxy
region$Status = ifelse(region$hg19POS %in% riskLoci$hg19POS1, "Index", "Proxy")
## What is the index snp for each row
region$IndexSNP = NA
region$IndexSNP_hg19POS = NA
for (i in 1:nrow(region)) {
	proxInd = which(riskLoci$hg19POS2 == region$hg19POS[i]) ## row of proxy
	region$IndexSNP[i] = riskLoci$SNP1_Name[proxInd]
	region$IndexSNP_hg19POS[i] = riskLoci$hg19POS1[proxInd]
}
## was index snp checked in eqtl analysis at all
region$IndexSNP_indata = NA
for (i in 1:nrow(region)) {
	indexInd = which(riskLoci_full$hg19POS2 == region$IndexSNP_hg19POS[i]) ## row of proxy
	region$IndexSNP_indata[i] = riskLoci_full$SNP2_missing[indexInd]
}
## does index snp have any significant eqtl result
region$IndexSNP_fdrSig = "na"
for (i in 1:nrow(region)) {
	if (region$IndexSNP_indata[i] == "analyzed") {
		pos = region$IndexSNP_hg19POS[i]
		sig = which(region$hg19POS == pos)
		region$IndexSNP_fdrSig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
	}
}
## what is the most significant eqtl result for the index snp
region$IndexSNP_mostSigFeat = NA
region$IndexSNP_mostSigFeat_gene = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {
		pos = region$IndexSNP_hg19POS[i]
		tmp = region[which(region$hg19POS == pos),]
		tmp = tmp[order(tmp$pvalue, decreasing=FALSE),]
		region$IndexSNP_mostSigFeat[i] = as.character(tmp$gene[1])
		region$IndexSNP_mostSigFeat_gene[i] = as.character(tmp$Symbol[1])
	}
}
## what is the lead variant for each SNP -- index if sig, else highest LD proxy
region$leadVariant = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {  ## index snp is fdr significant
		region$leadVariant[i] = region$IndexSNP[i]
		
	} else {									## index snp is NOT fdr significant
		## find highest LD proxy snp that is sig
		pos = region$IndexSNP_hg19POS[i]
		tmp = riskLoci[which(riskLoci$hg19POS1 == pos),]
		t_ind = which(tmp$Distance==0 | tmp$proxy_FDRsig==FALSE)
		if (length(t_ind>0)) { tmp = tmp[-t_ind,] }
		tmp = tmp[order(tmp$R_squared, rev(tmp$Distance), decreasing=TRUE),]
		region$leadVariant[i] = as.character(tmp$SNP2_Name[1])
		}
}
region$leadVariant_indicator = (region$SNP == region$leadVariant)

amyg = region


write.csv(sacc, file="raggr_31_snps_sacc_eqtls_fdr01.csv")
write.csv(amyg, file="raggr_31_snps_amyg_eqtls_fdr01.csv")
write.csv(dlp, file="raggr_31_snps_dlpfc_eqtls_fdr01.csv")



























#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("mergedEqtl_output_sacc_raggr_4features.rda", verbose=TRUE)
sigEqtlSacc_sub = allEqtl[allEqtl$FDR < 0.01,]
load("mergedEqtl_output_amyg_raggr_4features.rda", verbose=TRUE)
sigEqtlAmy_sub = allEqtl[allEqtl$FDR < 0.01,]

load("mergedEqtl_output_dlpfc_raggr_4features.rda", verbose=TRUE)
sigEqtlDlpfc_sub = allEqtl[allEqtl$FDR < 0.01,]

################
## metrics

## total features
nrow(sigEqtlSacc_sub)  ## 38609
nrow(sigEqtlAmy_sub)   ## 22572
nrow(sigEqtlDlpfc_sub) ## 17539

## per feature
table(sigEqtlSacc_sub$Type)
# Exon 	Gene  Jxn   Tx
# 20014  4463  8641  5491
table(sigEqtlAmy_sub$Type)
# Exon 	 Gene  Jxn   Tx
# 10370  2638  6435  3129
table(sigEqtlDlpfc_sub$Type)
# Exon  Gene  Jxn   Tx
# 7303  2381  4272  3583

## unique ensemblIDs
tapply(sigEqtlSacc_sub$EnsemblGeneID, sigEqtlSacc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  234  175  177  160
tapply(sigEqtlAmy_sub$EnsemblGeneID, sigEqtlAmy_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  145  105  130  105
tapply(sigEqtlDlpfc_sub$EnsemblGeneID, sigEqtlDlpfc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  136  105  108  127


################
## make csv

##### amygdala and sACC #####
dlp = sigEqtlDlpfc_sub
amyg = sigEqtlAmy_sub
sacc = sigEqtlSacc_sub
amyg$EnsemblGeneID = ss(amyg$EnsemblGeneID, "\\.")
sacc$EnsemblGeneID = ss(sacc$EnsemblGeneID, "\\.")

## snpMap
load("../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
gMap = as.data.frame(rowRanges(rse_gene))[,c("seqnames","start","end","strand","Class")]
eMap = as.data.frame(rowRanges(rse_exon))[,c("seqnames","start","end","strand","Class")]
jMap = as.data.frame(rowRanges(rse_jxn))[,c("seqnames","start","end","strand","Class")]
txMap = as.data.frame(rowRanges(rse_tx))[,c("seqnames","start","end","strand","source")]
txMap$source = "InGen"
# rm(rse_gene, rse_exon, rse_jxn, rse_tx)
colnames(gMap) = colnames(eMap) = colnames(jMap) = colnames(txMap) = 
	c("feat_chr","feat_start","feat_end","strand","Class")
featMap = rbind(rbind(rbind(gMap, eMap),jMap),txMap)
featMap$Type = c(rep("Gene",nrow(gMap)),rep("Exon",nrow(eMap)),rep("Jxn",nrow(jMap)),rep("Tx",nrow(txMap)))

geneMap = as.data.frame(rowRanges(rse_gene))[,c("gencodeID","Symbol","ensemblID","gene_type")]

## put together
snpMap_temp = snpMap[amyg$snps,]
featMap_temp = featMap[amyg$gene,]
geneMap_temp = geneMap[match(amyg$EnsemblGeneID, geneMap$ensemblID),]
amyg2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),amyg)
amyg3 = amyg2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(amyg3, "raggr_suggestive881_snps_amyg_eqtls_fdr01.csv")

snpMap_temp = snpMap[sacc$snps,]
featMap_temp = featMap[sacc$gene,]
geneMap_temp = geneMap[match(sacc$EnsemblGeneID, geneMap$ensemblID),]
sacc2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),sacc)
sacc3 = sacc2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(sacc3, "raggr_suggestive881_snps_sacc_eqtls_fdr01.csv")

tmp = amyg3
tmp$Status = NA
tmp$IndexSNP = NA
tmp$IndexSNP_hg19POS = NA
tmp$IndexPGCsig = NA
tmp$IndexFDR = NA


##### DLPFC #####
dlp = sigEqtlDlpfc_sub
dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")

## load SNP data
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda")

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
gMap = as.data.frame(rowRanges(rse_gene))[,c("seqnames","start","end","strand","Class")]
eMap = as.data.frame(rowRanges(rse_exon))[,c("seqnames","start","end","strand","Class")]
jMap = as.data.frame(rowRanges(rse_jxn))[,c("seqnames","start","end","strand","Class")]
txMap = as.data.frame(rowRanges(rse_tx))[,c("seqnames","start","end","strand","source")]
txMap$source = "InGen"
rm(rse_gene, rse_exon, rse_jxn, rse_tx)
colnames(gMap) = colnames(eMap) = colnames(jMap) = colnames(txMap) = 
	c("feat_chr","feat_start","feat_end","strand","Class")
featMap = rbind(rbind(rbind(gMap, eMap),jMap),txMap)
featMap$Type = c(rep("Gene",nrow(gMap)),rep("Exon",nrow(eMap)),rep("Jxn",nrow(jMap)),rep("Tx",nrow(txMap)))

geneMap = as.data.frame(rowRanges(rse_gene))[,c("gencodeID","Symbol","ensemblID","gene_type")]


snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
dlp2 = dlp2[,-which(colnames(dlp2)=="gencodeTx")]

dlp3 = dlp2[,c(2,12,13,25,19,14:18,21:23,26:29)]
write.csv(dlp3, "raggr_suggestive881_snps_dlpfc_eqtls_fdr01.csv")







