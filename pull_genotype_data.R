########################

## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)

## load data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda")
pd = colData(rse_gene)

## get BrNums
BrNums = BrNumOriginal = pd$BrNum
BrUniqueOriginal = unique(pd$BrNum)

### temporarily mislabel BrNumbers to pull correct genotype info
# - swap DNA BrNum labels from Br2260 and Br2473 when pulling genotype data
BrNums[which(BrNums %in% c("Br2260","Br2473"))]
BrNums[which(BrNums %in% c("Br2260","Br2473"))] = c("Br2473","Br2260","Br2473","Br2260")
BrNums[which(BrNums %in% c("Br2260","Br2473"))]
# - swap DNA BrNum labels from Br2301 and Br2538 when pulling genotype data
BrNums[which(BrNums %in% c("Br2301","Br2538"))]
BrNums[which(BrNums %in% c("Br2301","Br2538"))] = c("Br2538","Br2301","Br2301")
BrNums[which(BrNums %in% c("Br2301","Br2538"))]
# - change DNA label from Br2385 to Br2533 when pulling genotypes
BrNums[which(BrNums %in% c("Br2385","Br2533"))]
BrNums[which(BrNums %in% c("Br2385","Br2533"))] = c("Br2385","Br2385")
BrNums[which(BrNums %in% c("Br2385","Br2533"))]
# - swap DNA BrNum labels from Br5434 and Br5435 when pulling genotype data
BrNums[which(BrNums %in% c("Br5434","Br5435"))]
BrNums[which(BrNums %in% c("Br5434","Br5435"))] = c("Br5435","Br5434")
BrNums[which(BrNums %in% c("Br5434","Br5435"))]

pd$BrNum = BrNums
BrUniqueSwapped = unique(BrNums)

############# 
# fam file ##

### read in fam
# fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10.fam",
    # as.is=TRUE)
fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10.fam",
    as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")
fam$BrNum = fam$FID
fam$BrNum[fam$FID == "Omni2pt5"] = fam$IID[fam$FID == "Omni2pt5"]
# fam[grep("^Br", fam$IID),1:2] = fam[grep("^Br", fam$IID),2:1]
table(pd$BrNum %in% fam$BrNum) # keep samples w/ genotypes
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M

famOut = fam[which(fam$BrNum %in% pd$BrNum),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
# bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10"
bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10"
newbfile = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6"



## extract
system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.1 --maf 0.05 --hwe 0.000001 --make-bed --out", newbfile))

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
snp = as.matrix(snp[,unique(BrNums)])

### update MAP
# load("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10_updatedMap.rda")
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
mds = mds[unique(BrNums),-(1:3)]
colnames(mds) = paste0("snpPC",1:ncol(mds))


##########################
## correct BrNumbers #####

## confirm order stayed the same throughout
identical(rownames(mds), colnames(snp))
identical(rownames(mds), BrUniqueSwapped)

## reset to original correct BrNums
rownames(mds) = colnames(snp) = BrUniqueOriginal
pd$BrNum = BrNumOriginal



#############
## save #####
save(mds, snp, snpMap, compress=TRUE,
	file = "genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
save(mds, compress=TRUE,
	file = "genotype_data/zandiHyde_bipolar_MDS_n511.rda")
	
	
