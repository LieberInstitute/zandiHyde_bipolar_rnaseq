# module load plink/1.90b6.6 ## not needed

## Adapted from: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/pull_genotype_data.R

## Copy the files:
# cp ../genotype_data/zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6.bed zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6_unswaped.bed
# cp ../genotype_data/zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6.bim zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6_unswaped.bim

library('SummarizedExperiment')
library('sessioninfo')


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

## Read in the FAM file that has the swaps
fam = read.table("../genotype_data/zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6.fam", as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")

## Un-swap
bnums <- data.frame(
    Fixed  = BrNumOriginal,
    Mislabeled = BrNums,
    stringsAsFactors = FALSE
)
bnums$diff  <- bnums$Fixed != bnums$Mislabeled
table(bnums$diff)
#
# FALSE  TRUE
#   500    11

subset(bnums, diff)
#      Fixed Mislabeled diff
# 153 Br2533     Br2385 TRUE
# 167 Br2260     Br2473 TRUE
# 173 Br2473     Br2260 TRUE
# 176 Br5434     Br5435 TRUE
# 181 Br5435     Br5434 TRUE
# 196 Br2301     Br2538 TRUE
# 362 Br2260     Br2473 TRUE
# 367 Br2473     Br2260 TRUE
# 377 Br2533     Br2385 TRUE
# 381 Br2538     Br2301 TRUE
# 414 Br2538     Br2301 TRUE

table(fam$FID %in% bnums$Mislabeled)
#
# TRUE
#  295

bnums_diff <- subset(bnums, diff)
table(fam$FID %in% bnums_diff$Mislabeled)
# FALSE  TRUE
#   288     7
length(unique(bnums_diff$Mislabeled))
# [1] 7

tofix <- fam$FID %in% bnums_diff$Mislabeled
m <- match(fam$FID[tofix], bnums_diff$Mislabeled)

fam_fixed <- fam
fam_fixed$FID[tofix] <- bnums_diff$Fixed[m]
table(fam_fixed$FID == fam$FID)
# FALSE  TRUE
#     7   288

write.table(fam_fixed, file = 'zandiHyde_bipolar_Genotypes_n511_maf005_geno10_hwe1e6_unswaped.fam', col.names = FALSE, sep = ' ', quote = FALSE, row.names = FALSE)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-12-05
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.0    2019-10-30 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  digest                 0.6.22    2019-10-21 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggplot2                3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.0    2019-10-29 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  rlang                  0.4.1     2019-10-24 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  S4Vectors            * 0.24.0    2019-10-29 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.0    2019-10-29 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.10      2019-10-01 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library
