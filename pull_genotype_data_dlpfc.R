########################

######################## get 10,777 snps used in other regions
## load SNP data
load("genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("eqtl_raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snpMap_10777 = snpMap

########################
## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)

## load samples
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
BrainNums = colData(rse_gene)$BrNum
rm(rse_gene)

############# 
# fam file ##

### read in fam
fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10.fam",
    as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")
fam$BrNum = fam$FID
fam$BrNum[fam$FID == "Omni2pt5"] = fam$IID[fam$FID == "Omni2pt5"]
table(BrainNums %in% fam$BrNum) # keep samples w/ genotypes
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M

famOut = fam[which(fam$BrNum %in% BrainNums),]
write.table(famOut[,1:2], "samples_to_extract167.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(rownames(snpMap_10777), "snps_to_extract10777.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10" ## keep
newbfile = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777"   ## change

# https://www.staff.ncl.ac.uk/heather.cordell/mres2012plink.html

## extract
system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	"--keep samples_to_extract167.txt --extract snps_to_extract10777.txt --make-bed --out ", newbfile))

# ## extract
# system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	# "--keep samples_to_extract167.txt --geno 0.01 --maf 0.005 --hwe 0.00000001 --make-bed --out", newbfile))

# ## independent and cluster
system(paste("/users/ajaffe/bin/plink --bfile", newbfile, "--indep 100 10 1.25 --out", newbfile))

## MDS components	
system(paste0("/users/ajaffe/bin/plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))

# ## A transpose
system(paste("/users/ajaffe/bin/plink --bfile", newbfile,
	"--recode A-transpose --out", newbfile))
	
################
## read in #####

## read in genotypes
genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")

snp = as.data.frame(genotypes[,-(1:6)])
colnames(snp) = ifelse(grepl("^Br", ss(colnames(snp), "_")),
			ss(colnames(snp), "_"), ss(colnames(snp), "_",2))
snp = as.matrix(snp[,unique(BrainNums)])

### update MAP
load("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10_updatedMap.rda")
snpMap = snpMap[match(genotypes$SNP, snpMap$SNP),]
rownames(snp) = rownames(snpMap) =  snpMap$SNP

## check directionality
swapIndex = which(snpMap$ALT == genotypes$COUNTED)
snp[swapIndex,] = 2-snp[swapIndex,]

#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE)
rownames(mds) = ifelse(grepl("^Br", mds$FID),
			mds$FID, mds$IID)
mds = mds[unique(BrainNums),-(1:3)]
colnames(mds) = paste0("snpPC",1:ncol(mds))


#############
## save #####
save(mds, snp, snpMap, compress=TRUE,
	file = "genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
	
	
