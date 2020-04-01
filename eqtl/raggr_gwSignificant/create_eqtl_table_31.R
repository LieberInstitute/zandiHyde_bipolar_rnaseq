##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


snpspos$downstream = snpspos$pos-5e5
snpspos$upstream = snpspos$pos+5e5
snpspos$downstream[snpspos$downstream < 0] = 0


which(posGene$seqnames == snpspos$chr & 											## same chr
		( (posGene$start>snpspos$downstream & posGene$start<snpspos$upstream) | 	# starts inside window OR
		  (posGene$end>snpspos$downstream & posGene$end<snpspos$upstream) |			# ends inside window OR
		  (posGene$start<snpspos$downstream & posGene$end>snpspos$upstream) ) )		# covers whole window

#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("../raggr/mergedEqtl_output_sacc_raggr_4features.rda", verbose=TRUE)
sigEqtlSacc = allEqtl
load("../raggr/mergedEqtl_output_amyg_raggr_4features.rda", verbose=TRUE)
sigEqtlAmy = allEqtl
load("../raggr/mergedEqtl_output_dlpfc_raggr_4features.rda", verbose=TRUE)
sigEqtlDlpfc = allEqtl

################
## load SNPs
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## add coordinate to eqtl results
snpMap2 = snpMap[(snpMap$SNP %in% sigEqtlSacc$snps),]
posInd = match(sigEqtlSacc$snps, snpMap2$SNP)
sigEqtlSacc$hg19POS = snpMap2$pos_hg19[posInd]

snpMap2 = snpMap[(snpMap$SNP %in% sigEqtlAmy$snps),]
posInd = match(sigEqtlAmy$snps, snpMap2$SNP)
sigEqtlAmy$hg19POS = snpMap2$pos_hg19[posInd]

snpMap2 = snpMap[(snpMap$SNP %in% sigEqtlDlpfc$snps),]
posInd = match(sigEqtlDlpfc$snps, snpMap2$SNP)
sigEqtlDlpfc$hg19POS = snpMap2$pos_hg19[posInd]


################
## subset to 31 significant index SNPs
indexLoci = read.csv("../../PGC_risk_loci.csv", stringsAsFactors=FALSE) ## 881
## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## subset riskLoci table to rows with g.w. significant index SNP 
riskLoci$genomewideSig1 = NA
for (i in 1:nrow(riskLoci)) {
	indexInd = which(indexLoci$hg19POS == riskLoci$hg19POS1[i])
	riskLoci$genomewideSig1[i] = indexLoci$genomewide[indexInd]
}
riskLociGW = riskLoci[which(riskLoci$genomewideSig1==TRUE),]

# ## subset results
# sigEqtlSacc2 = sigEqtlSacc[which(sigEqtlSacc$snps %in% riskLociGW$SNP2_Name),]
# sigEqtlAmy2 = sigEqtlAmy[which(sigEqtlAmy$snps %in% riskLociGW$SNP2_Name),]
# sigEqtlDlpfc2 = sigEqtlDlpfc[which(sigEqtlDlpfc$snps %in% riskLociGW$SNP2_Name),]
## subset results
sigEqtlSacc = sigEqtlSacc[which(sigEqtlSacc$hg19POS %in% riskLociGW$hg19POS2),]
sigEqtlAmy = sigEqtlAmy[which(sigEqtlAmy$hg19POS %in% riskLociGW$hg19POS2),]
sigEqtlDlpfc = sigEqtlDlpfc[which(sigEqtlDlpfc$hg19POS %in% riskLociGW$hg19POS2),]


## recalculate FDR
sigEqtlSacc$FDR = p.adjust(sigEqtlSacc$pvalue, method="fdr")
sigEqtlAmy$FDR = p.adjust(sigEqtlAmy$pvalue, method="fdr")
sigEqtlDlpfc$FDR = p.adjust(sigEqtlDlpfc$pvalue, method="fdr")

write.csv(sigEqtlAmy, file="all_raggr_31_snps_amyg_eqtls.csv")
write.csv(sigEqtlSacc, file="all_raggr_31_snps_sacc_eqtls.csv")
write.csv(sigEqtlDlpfc, file="all_raggr_31_snps_dlpfc_eqtls.csv")


######################################
## start from beginning of original script

sigEqtlSacc_sub = sigEqtlSacc[sigEqtlSacc$FDR < 0.01,]
sigEqtlAmy_sub = sigEqtlAmy[sigEqtlAmy$FDR < 0.01,]
sigEqtlDlpfc_sub = sigEqtlDlpfc[sigEqtlDlpfc$FDR < 0.01,]


################
## make csv

##### amygdala and sACC #####
dlp = sigEqtlDlpfc_sub
amyg = sigEqtlAmy_sub
sacc = sigEqtlSacc_sub
amyg$EnsemblGeneID = ss(amyg$EnsemblGeneID, "\\.")
sacc$EnsemblGeneID = ss(sacc$EnsemblGeneID, "\\.")

## snpMap
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load("../../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../../data/zandiHypde_bipolar_rseGene_n511.rda")
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
write.csv(amyg3, "raggr_genomewide31_snps_amyg_eqtls_fdr01.csv")

snpMap_temp = snpMap[sacc$snps,]
featMap_temp = featMap[sacc$gene,]
geneMap_temp = geneMap[match(sacc$EnsemblGeneID, geneMap$ensemblID),]
sacc2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),sacc)
sacc3 = sacc2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(sacc3, "raggr_genomewide31_snps_sacc_eqtls_fdr01.csv")

snpMap_amyg = snpMap
##### DLPFC #####
dlp = sigEqtlDlpfc_sub
dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")

## load SNP data
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda", verbose=TRUE)
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap2 = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
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


snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
dlp2 = dlp2[,-which(colnames(dlp2)=="gencodeTx")]

dlp3 = dlp2[,c(2,12,13,14,26,20,15:19,22:24,27:30)]

## fill in NA snps
missInd = which(is.na(dlp3$SNP))
missSnp = dlp$snps[missInd]
snpMap_miss = snpMap_amyg[missSnp,]
count = 1
for (i in missInd) {  
	rownames(dlp3)[i] = rownames(snpMap_miss)[count]
	dlp3[i,1:4] = c(snpMap_miss[count,1:4])
	count = count+1
}

write.csv(dlp3, "raggr_genomewide31_snps_dlpfc_eqtls_fdr01.csv")

