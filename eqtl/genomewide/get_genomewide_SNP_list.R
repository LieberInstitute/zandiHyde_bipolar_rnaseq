########################

######################## get overlapping snps
## load SNP data

## Amyg / sACC snps
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
snpMap_amyg = snpMap

## DLPFC snps
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS
rm(snp, mds)


ind = which(snpMap$SNP %in% snpMap_amyg$SNP)
snpMap = snpMap[ind,]
snpMap_amyg = snpMap_amyg[rownames(snpMap2),]

# ## check
# identical(snpMap$POS, snpMap_amyg$POS)
# identical(snpMap$newRef, snpMap_amyg$newRef)
# identical(snpMap$rsNumGuess, snpMap_amyg$rsNumGuess)

snpMapKeep = snpMap
save(snpMapKeep, file="rdas/overlappingSNPs.rda")

