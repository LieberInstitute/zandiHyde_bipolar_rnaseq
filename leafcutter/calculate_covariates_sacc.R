####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
print("Loaded SACC dataframe")
### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"
print("Created Other dx bipolar dataframe")
## load SNP data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
print("Loaded SNP data")
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
snp = snp[-snpInd,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 
print("Filtered the brain region")
######################
# statistical model ##
######################
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))

mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]
print("Generated statistical model")
######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")
print("Created SNP object")
#######################
####### do PCA ########
#######################

# load the splice ratio file
splice <- read.table("/users/schadinh/lieber/qqnorm/sacc_qqnorm_20190220.txt.gz",sep='\t',comment.char='$',header=T)
# drop meta data columns
splice <- splice[,5:dim(splice)[2]]
# run PCA
sPCA <- prcomp(t(splice))
numPCs <- num.sv(splice, mod)
sPCs <- sPCA$x[,1:numPCs]
covsSplice <- t(cbind(mod[,-1],sPCs))
write.table(covsSplice,file="covariates_sACC.txt",sep="\t",quote=F,col.names=T,row.names=F)
