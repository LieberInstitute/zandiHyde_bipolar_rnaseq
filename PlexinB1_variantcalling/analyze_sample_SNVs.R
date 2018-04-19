##
library(jaffelab)
library(VariantAnnotation)

### pd file / sample IDs
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata")

### make "Other" dx bipolar
pd = colData(rse_gene)[,1:13]
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"
pd$PrimaryDx = droplevels(pd$PrimaryDx)

## confirm pd order matches vcf order
samples = read.table("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/samples.manifest", stringsAsFactors=FALSE)
identical(samples$V5, as.character(pd$SAMPLE_ID))


##########################
### Load in snpMap
snpMap = rtracklayer::import("plexinb1_hg38.bed")

##########################
# read in merged VCF file
mergedVcfFile = 'Genotypes/mergedVariants.vcf.gz'
vcf = readVcf(mergedVcfFile,'GRCh38.p7')
info(vcf)$RS = mcols(snpMap)$name[match(rowRanges(vcf),snpMap)]
colnames(vcf) = paste0(pd$PrimaryDx, "_", pd$RNum)


## sample indexes of 511 samples used in analyses
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda")
id511 = colData(rse_gene)
keepInd = which(pd$SAMPLE_ID %in% id511$SAMPLE_ID)
pd = pd[keepInd,]


########################################
# plot snp correlation of all samples
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snps = snps[,keepInd]

rownames(snps)[grep("48422302",rownames(snps))] = "chr3:48422302_C/T_rs9883927"
rownames(snps)[grep("48422390",rownames(snps))] = "chr3:48422390_T/G_POSofINTEREST"

pd2 = cbind(pd, t(snps))
for (i in 14:33) { print(table(pd2[,i], pd2$PrimaryDx)) }



########################################
# unique brains

keepInd = which(!duplicated(pd$BrNum))
pd = pd[keepInd,]
snps = snps[,keepInd]

snpList = list()
pd2 = cbind(pd, t(snps))
for (i in 14:33) { 
	pd2[,i] = factor(pd2[,i], levels=c("0","1","2"))
	snpList[[i-13]] = table(pd2[,i], pd2$PrimaryDx)
	names(snpList)[i-13] = colnames(pd2)[i]
	}
	
pd3 = pd2[which(pd2$"chr3.48422390_T.G_POSofINTEREST" != 0),]
pd3[order(pd3$PrimaryDx,pd3$BrNum),c(1,3,6:8,11,23)]













