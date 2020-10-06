library(data.table)

setDTthreads(1)

load("/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38.Rdata")

snpMap <- as.data.table(snpMap)

hg19_gwas <- fread("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/daner_PGC_BIP32b_mds7a_0416a")

snpMap$hg19_key <- paste0(snpMap$CHR, "_", snpMap$POS)

hg19_gwas$hg19_key <- paste0(hg19_gwas$CHR, "_", hg19_gwas$BP)

hg38_gwas <- merge(hg19_gwas, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

hg38_gwas <- hg38_gwas[,-c("hg19_key", "CHR", "BP")]

col_order <- c(
    "chr_hg38", "SNP", "pos_hg38", "A1", "A2", "FRQ_A_20352",
    "FRQ_U_31358", "INFO", "OR", "SE", "P", "ngt", "Direction",
    "HetISqt", "HetDf", "HetPVa", "Nca", "Nco", "Neff")

hg38_gwas <- hg38_gwas[, ..col_order]

names(hg38_gwas)[1] <- "CHR"

names(hg38_gwas)[3] <- "BP"

write.table(hg38_gwas, file = "PGC_BIP_hg38.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
