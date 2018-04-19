library(jaffelab)
library(VariantAnnotation)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
col.pal = brewer.pal(9,"Blues")

##########################
### pd file / sample IDs for n = 540
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata")
pd = colData(rse_gene)[,1:11]
pd$samplelabel = paste0(pd$BrNum,"_",pd$Brain.Region,"_",pd$RNum)

##########################
### Load in snpMap
snpMap = rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")

##########################
# read in merged VCF file
mergedVcfFile = '/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/Genotypes/mergedVariants.vcf.gz'
vcf = readVcf(mergedVcfFile,'GRCh38.p2')
info(vcf)$RS = mcols(snpMap)$name[match(rowRanges(vcf),snpMap)]
all(colnames(vcf) == colData(rse_gene)$bamFile)  ## check order
colnames(vcf) = pd$samplelabel


######################
# subset to high-depth
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
          info(vcf)$VDB >0.1,]


########################################
# plot snp correlation of all samples
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snpCor = cor(snps, use="pairwise.complete.obs")
rownames(snpCor) = colnames(snpCor) = pd$samplelabel


# # plot
# pdf("pheatmap_n540.pdf",h=40,w=40)
# pheatmap(snpCor, 
		# cluster_rows=T, 
		# cluster_cols=T,
		# color=col.pal)
# dev.off()



### Samples with high correlation, unmatching BrNums
l = list()
listInd = 1

for (i in 1:540) {
	matchInd = which(snpCor[i,]>0.6)
	if (length(matchInd) > 1) {
	subSnpCor = snpCor[matchInd,matchInd]
	ids = ss(colnames(subSnpCor),"_")
		if (length(unique(ids))>1) { 
			l[[listInd]] = subSnpCor
			listInd = listInd+1
		}
	}
}
length(unique(l))
unique(l)



### Samples with matching BrNums, low correlation
l = list()
listInd = 1
br = as.character(unique(pd$BrNum))

for (i in 1:length(br)) {
	matchInd = which(pd$BrNum == br[i])
	if (length(matchInd) > 1) {
	subSnpCor = snpCor[matchInd,matchInd]
	if (sum(subSnpCor < .6) > 0) {
		l[[listInd]] = subSnpCor
		listInd = listInd+1
		}
	}
}
length(l)
l
                       # Br5266_Amygdala_R13922 Br5266_sACC_R14138
# Br5266_Amygdala_R13922              1.0000000          0.3444686
# Br5266_sACC_R14138                  0.3444686          1.0000000

# [[2]]
                       # Br1697_Amygdala_R13976 Br1697_sACC_R14200
# Br1697_Amygdala_R13976              1.0000000          0.2915438
# Br1697_sACC_R14200                  0.2915438          1.0000000

# [[4]]
                       # Br5205_Amygdala_R14019 Br5205_sACC_R14244
# Br5205_Amygdala_R14019              1.0000000          0.3873371
# Br5205_sACC_R14244                  0.3873371          1.0000000

# [[6]]
                       # Br5168_Amygdala_R14059 Br5168_sACC_R14278
# Br5168_Amygdala_R14059              1.0000000          0.1580363
# Br5168_sACC_R14278                  0.1580363          1.0000000

# [[7]]
                       # Br5434_Amygdala_R14080 Br5434_sACC_R14299
# Br5434_Amygdala_R14080              1.0000000          0.3458916
# Br5434_sACC_R14299                  0.3458916          1.0000000



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Check c("Br5266","Br1697","Br5205","Br5168","Br5434") against older RNA data

Br = c("Br5266","Br1697","Br5205","Br5168","Br5434")
RNums = as.character(pd$RNum[pd$BrNum %in% Br])

##########################
# read in merged VCF file
mergedVcfFile = '/dcl01/lieber/ajaffe/lab/brain_swap/allMergedVariants.vcf.gz'
vcf = readVcf(mergedVcfFile,'GRCh38.p2')
info(vcf)$RS = mcols(snpMap)$name[match(rowRanges(vcf),snpMap)]
colnames(vcf) = ss(basename(colnames(vcf)),"\\.")

######################
# subset to high-depth
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
          info(vcf)$VDB >0.1,]

#################
# snp correlation
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snpCor = cor(snps, use="pairwise.complete.obs")
snpCorInt = snpCor[,grep("13976|14278",colnames(snpCor))]

matches = which(snpCorInt>.6, arr.ind=T)
matches = snpCorInt[unique(matches[,1]),unique(matches[,2])]
matches = round(matches,3)
colnames(matches) = ss(colnames(matches),"_")
matches

pd[pd$RNum%in%colnames(matches),]

## check BrNums of matching RNums 
load("/dcl01/lieber/ajaffe/lab/brain_swap/genotype_match_matrix_allBrains.rda")

snpCorBr1697 = snpCor2[,c(3443,3652,457,746,3143,3144)]
matches = which(snpCorBr > .6, arr.ind=TRUE)
snpCorBr[unique(matches[,1]),]

snpCorBr5168 = snpCor2[,c(3519,3727,2664,2665)]
matches = which(snpCorBr5168 > .7, arr.ind=TRUE)
snpCorBr5168[unique(matches[,1]),]

snpCorBr5205 = snpCor2[,c(561,985,1336,2672,2673,3482,3694)]
matches = which(snpCorBr5205 > .7, arr.ind=TRUE)
snpCorBr5205[unique(matches[,1]),]

snpCorOther = snpCor2[,c(3390,3591,3539,3747)]
matches = which(snpCorOther > .7, arr.ind=TRUE)
snpCorOther[unique(matches[,1]),]



snpCorBr = snpCor2[,grep("13922|14138|14080|14299",colnames(snpCor2))]
matches = which(snpCorBr > .6, arr.ind=TRUE)
snpCorBr[unique(matches[,1]),]


