####
### libraries
library(SummarizedExperiment)
library(jaffelab)

## risk loci from schiz PGC paper
szLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", 
					stringsAsFactors=FALSE)
szLoci = szLoci[,c(4:9, 25:28,24)]
szLoci$hg19POS = paste0(szLoci$snp_chr, ":", szLoci$snp_pos_hg19) 

## risk loci from bipolar PGC paper
bpLoci = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/rAggr_results_881.csv", 
					stringsAsFactors=FALSE)
colnames(bpLoci) = gsub("\\.", "_", colnames(bpLoci))
bpLoci$hg19POS1 = paste0(bpLoci$SNP1_Chr, ":", bpLoci$SNP1_Pos) 
bpLoci$hg19POS2 = paste0(bpLoci$SNP2_Chr, ":", bpLoci$SNP2_Pos) 

indexLoci = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/PGC_risk_loci.csv", 
					stringsAsFactors=FALSE) ## 881
bpLoci$shortSNP1 = ss(bpLoci$SNP1_Name, ":")
bpLoci$shortSNP2 = ss(bpLoci$SNP2_Name, ":")
bpLoci$Status1 = ifelse(bpLoci$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
bpLoci$Status2 = ifelse(bpLoci$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")

## test by both position and name
bpLoci$Scz2 = bpLoci$hg19POS2 %in% szLoci$hg19POS
bpLoci$Scz22 = bpLoci$shortSNP2 %in% szLoci$Index_SNP_dbSNP_b141_

indexLoci$Scz = indexLoci$hg19POS %in% szLoci$hg19POS
indexLoci$Scz1 = indexLoci$SNP %in% szLoci$Index_SNP_dbSNP_b141_

b = bpLoci[which(bpLoci$Scz2 == TRUE | bpLoci$Scz22 == TRUE),]
i = indexLoci[which(indexLoci$Scz== TRUE | indexLoci$Scz1 == TRUE),]


loci = unique(b$shortSNP2)
# [1] "rs140505938" "rs77011057"  "rs13169274"  "rs12704290"  "rs7893279"   "rs9545047"   "rs35604463"  "rs12908161"

szLoci[which(szLoci$Index_SNP_dbSNP_b141_ %in% loci),]
bpLoci[which(bpLoci$shortSNP2 %in% loci),]

# combine, make csv


scz = read.csv("schiz_pgc_loci_overlap.csv", stringsAsFactors=FALSE)

amyg = read.csv("conditional/raggr_881_snps_amyg_eqtls_fdr01.csv", row.names=1)
sacc = read.csv("conditional/raggr_881_snps_sacc_eqtls_fdr01.csv", row.names=1)
dlp = read.csv("conditional/raggr_881_snps_dlpfc_eqtls_fdr01.csv", row.names=1)

a = amyg[amyg$hg19POS %in% scz$hg19POS,]
s = sacc[sacc$hg19POS %in% scz$hg19POS,]
d = dlp[dlp$hg19POS %in% scz$hg19POS,]

library(sheetr)

# Create a list of dataframes
dataframes = list()
dataframes[["DLPFC"]] = d
dataframes[["Amygdala"]] = a
dataframes[["sACC"]] = s

# Write the list of dataframes to CSV file
write_dataframes_to_csv(dataframes, "schiz_pgc_loci_bipolar_results.csv")














	  