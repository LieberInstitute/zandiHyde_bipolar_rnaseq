library(data.table)
library(dplyr)

setDTthreads(1)

# Load snpMap
load(
    "/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38.Rdata"
)

snpMap <- as.data.table(snpMap)

# Read in Angela's GWAS
hg19_gwas <-
    fread(
        "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/daner_PGC_BIP32b_mds7a_0416a"
    )

# Make key
snpMap$hg19_key <- paste0(snpMap$CHR, "_", snpMap$POS)

hg19_gwas$hg19_key <- paste0(hg19_gwas$CHR, "_", hg19_gwas$BP)

# Merge in the hg38 coordinates based on hg19 coordinates
hg38_gwas <-
    merge(hg19_gwas, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

# Drop hg19 coords from gwas
hg38_gwas <- hg38_gwas[, -c("hg19_key", "CHR", "BP")]

# reorder columns
col_order <- c(
    "chr_hg38",
    "SNP",
    "pos_hg38",
    "A1",
    "A2",
    "FRQ_A_20352",
    "FRQ_U_31358",
    "INFO",
    "OR",
    "SE",
    "P",
    "ngt",
    "Direction",
    "HetISqt",
    "HetDf",
    "HetPVa",
    "Nca",
    "Nco",
    "Neff"
)

hg38_gwas <- hg38_gwas[, ..col_order]

names(hg38_gwas)[1] <- "CHR"

names(hg38_gwas)[3] <- "BP"

# read bim file with unique rsIDs
uniq_bim <-
    fread(
        "filter_data/unique_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_amygdala_gene_uniqueSNPs.bim"
    )

colnames(uniq_bim) <- c("CHR", "SNP", "dummy", "BP", "A1", "A2")

# check for overlap between gwas and unique bim file
uniq_bim$new_test <- paste0(uniq_bim$CHR, "_", uniq_bim$BP)

hg38_gwas$new_test <-
    paste0(gsub("chr", "" , hg38_gwas$CHR), "_", hg38_gwas$BP)

table(uniq_bim$new_test %in% hg38_gwas$new_test)

# > table(uniq_bim$new_test %in% hg38_gwas$new_test)
#
#   FALSE    TRUE
# 3101196 7885983

# merge in unique rsIDs
hg38_gwas <-
    merge(hg38_gwas[, -"SNP"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas$new_test <- NULL

col_order <- c(
    "CHR",
    "SNP",
    "BP",
    "A1",
    "A2",
    "FRQ_A_20352",
    "FRQ_U_31358",
    "INFO",
    "OR",
    "SE",
    "P",
    "ngt",
    "Direction",
    "HetISqt",
    "HetDf",
    "HetPVa",
    "Nca",
    "Nco",
    "Neff"
)

hg38_gwas <- hg38_gwas[, ..col_order]

hg38_gwas$effect <- log(hg38_gwas$OR)
hg38_gwas$Z <- hg38_gwas$effect / hg38_gwas$SE

pdf(file = "PGC_BIP_hg38_hist.png.pdf", useDingbats = FALSE)

hist(hg38_gwas$effect, color = "gold")
hist(hg38_gwas$Z, color = "darkorange")

dev.off()

hg38_gwas_clean <-
    hg38_gwas[, c("SNP", "A1", "A2", "Neff", "P", "Z")]
colnames(hg38_gwas_clean)[4] <- "N"

write.table(
    hg38_gwas_clean,
    file = "PGC_BIP_hg38_clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
