##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)

################
## load SNP data
load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## risk loci from PGC paper
indexLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors=FALSE) ## 179
indexLoci$hg19POS = paste0(indexLoci$Chromosome, ":", indexLoci$snp_pos_hg19)
indexIndex = which(snpMap$pos_hg19 %in% indexLoci$hg19POS)	# keep 140

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from raggr list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 9798 snps in snpMap
snpMap = snpMap[keepIndex,]
keepIndex = which(riskLoci$hg19POS2 %in% snpMap$pos_hg19)	# keep 9798 snps in riskLoci
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
# 179 indexSNPs
# 10,981 rAggr riskLoci (including index SNPs)  		#  length(unique(riskLoci_full$SNP2_Name))
# 140 indexSNPs ~IN~ our SNPs used in eQTLS				#  length(which(snpMap$pos_hg19 %in% riskLoci$hg19POS1))
# 164 indexSNPs ~OF PROXIES~ in our SNPs used in eQTLS	#  length(unique(riskLoci$hg19POS1))
# 9,798 rAggr riskLoci in our SNPs used in eQTLS 		#  length(unique(riskLoci$hg19POS2))


################
## load table
hippo = read.csv("raggr_suggestive881_snps_hippo_eqtls_fdr05.csv", row.names=1)
dg = read.csv("raggr_suggestive881_snps_dg_eqtls_fdr05.csv", row.names=1)
inter = read.csv("raggr_suggestive881_snps_inter_eqtls_fdr05.csv", row.names=1)


## unique SNPs
length(unique(dg$hg19POS))		# 3286
length(unique(hippo$hg19POS))	# 2898
length(unique(inter$hg19POS))	# 563
v = venn.diagram(list(HIPPO = hippo$hg19POS, DG = dg$hg19POS),
	fill = c("lightgoldenrod", "coral"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_SNP.pdf", useDingbats = FALSE)
    grid.draw(v)
dev.off()

## unique index SNPs
length(unique(dg$hg19POS[dg$hg19POS %in% riskLoci$hg19POS1]))			# 57
length(unique(hippo$hg19POS[hippo$hg19POS %in% riskLoci$hg19POS1]))		# 45
length(unique(inter$hg19POS[inter$hg19POS %in% riskLoci$hg19POS1]))		# 9
v = venn.diagram(list(HIPPO = unique(hippo$hg19POS[hippo$hg19POS %in% riskLoci$hg19POS1]), 
				DG = unique(dg$hg19POS[dg$hg19POS %in% riskLoci$hg19POS1])), 
	fill = c("lightgoldenrod", "coral"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_index.pdf", useDingbats = FALSE)
    grid.draw(v)
dev.off()
	
## unique features
tapply(dg$gene, dg$Type, function(x) length(unique(x)))
tapply(hippo$gene, hippo$Type, function(x) length(unique(x)))
tapply(inter$gene, inter$Type, function(x) length(unique(x)))
v = venn.diagram(list(HIPPO = hippo$gene, DG = dg$gene), 
	fill = c("lightgoldenrod", "coral"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_features.pdf", useDingbats = FALSE)
    grid.draw(v)
dev.off()

## SNP-feature pairs
nrow(dg)  		## 25410
nrow(hippo)		## 10914
nrow(inter)		## 984
table(dg$Type)
table(hippo$Type)
table(inter$Type)
v = venn.diagram(list(HIPPO = paste0(hippo$SNP, hippo$gene),  DG = paste0(dg$SNP, dg$gene)), 
	fill = c("lightgoldenrod", "coral"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_snpfeatpairs.pdf", useDingbats = FALSE)
    grid.draw(v)
dev.off()

## Unique symbols in SNP-feature pairs
length(unique(dg$Symbol))		# 233
length(unique(hippo$Symbol))	# 148
length(unique(inter$Symbol))	# 19
tapply(dg$Symbol, dg$Type, function(x) length(unique(x)))
tapply(hippo$Symbol, hippo$Type, function(x) length(unique(x)))
tapply(inter$Symbol, inter$Type, function(x) length(unique(x)))
v = venn.diagram(list(HIPPO = hippo$Symbol[-is.na(hippo$Symbol)], DG = dg$Symbol[-is.na(dg$Symbol)]), 
	fill = c("lightgoldenrod", "coral"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_symbol.pdf", useDingbats = FALSE)
    grid.draw(v)
dev.off()
	
	
	
	
################################################################
###### Index SNP info ##########
################################################################

region = inter

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
proxInd = match(region$hg19POS, riskLoci$hg19POS2)
region$IndexSNP = riskLoci$SNP1_Name[proxInd]
region$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]

## was index snp checked in eqtl analysis at all
indexInd = match(region$IndexSNP_hg19POS, riskLoci_full$hg19POS2) ## row of proxy
region$IndexSNP_indata = riskLoci_full$SNP2_missing[indexInd]

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
## Is the SNP the lead variant? (Check by position)
leadVarInd = match(region$leadVariant, riskLoci$SNP2_Name)
leadVarPos = riskLoci$hg19POS2[leadVarInd]
region$leadVariant_indicator = (region$hg19POS == leadVarPos)


inter = region


write.csv(hippo, file="raggr_179_snps_hippo_eqtls_fdr05.csv")
write.csv(dg, file="raggr_179_snps_dg_eqtls_fdr05.csv")
write.csv(inter, file="raggr_179_snps_inter_eqtls_fdr05.csv")





## Index SNP numbers:
# HIPPO
# Index SNPs in our data with FDR result:
# 45		# length(unique(hippo[which(hippo$Status == "Index"),]$IndexSNP))
# Index SNPs from proxies with FDR result:
# 76		# length(unique(hippo$IndexSNP)), length(unique(hippo$leadVariant)

# DG
# Index SNPs in our data with FDR result
# 57		# length(unique(dg[which(dg$Status == "Index"),]$IndexSNP))
# Index SNPs from proxies with FDR result
# 103		# length(unique(dg$IndexSNP))

# Interaction
# Index SNPs in our data with FDR result:
# 9			# length(unique(inter[which(inter$Status == "Index"),]$IndexSNP))
# Index SNPs from proxies with FDR result:
# 21		# length(unique(inter$IndexSNP))










