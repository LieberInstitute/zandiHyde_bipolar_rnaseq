## Based on:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/read_twas.R
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas.R
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas_psycm.R


library("readr")
library("purrr")
library("dplyr")
library("stringr")
library("sessioninfo")
library("SummarizedExperiment")
library("ggplot2")
library("gplots")
library("VennDiagram")
library("RColorBrewer")
library("readr")
library("here")

## Main options to work through
regions <- c("amygdala", "sacc")
types <- c("pgc")
features <- "gene"
file_types <- c("all", "included", "dropped")

## Function for locating the files
locate_files <- function(ftype, path, type) {

    ## Determine the file pattern to use
    fpatt <- case_when(
        ftype == "all" ~ "",
        ftype == "included" ~ "\\.analysis\\.joint_included",
        ftype == "dropped" ~ "\\.analysis\\.joint_dropped"
    )
    patt <- paste0(
        type,
        "\\.[[:digit:]]*",
        fpatt,
        "\\.dat$"
    )

    ## Find the files
    dat_files <- dir(
        path = path,
        pattern = patt,
        full.names = TRUE
    )

    ## Extract the chromosome
    names(dat_files) <- str_extract(
        basename(dat_files),
        "[:digit:]+"
    )

    ## Done
    return(dat_files)
}


## Now read in the data for all file types
## Note that each file type has different numbers of columns
## when why I'm keeping them as different elements of the
## 'twas' list
twas <- map(file_types, function(ftype) {

    ## data.frame with all the combinations of arguments
    arg_grid <- expand.grid(
        region = regions,
        type = types,
        feature = features,
        stringsAsFactors = FALSE
    )


    pmap_dfr(arg_grid, function(region, type, feature) {

        ## Construct the path to the files
        path <- file.path(
            "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas",
            paste0(region,"_",feature),
            type
        )

        ## Now locate the files
        dat_files <- locate_files(
            ftype = ftype,
            path = path,
            type = type
        )

        ## Next read the files and add the chromosome info
        result <- map2_dfr(
            dat_files,
            names(dat_files),
            function(f, chr) {
                res <- suppressMessages(read_tsv(f))
                res$chr <- chr
                return(res)
            }
        )

        ## Next add the region, feature and type information
        result$region <- region
        result$feature <- feature
        result$type <- type

        ## Done
        return(result)
    })
})
names(twas) <- file_types

## Explore the resulting dimensions
map_dfr(twas, dim)
# # A tibble: 2 x 3
#     all included dropped
#   <int>    <int>   <int>
# 1 31609       79   16136
# 2    24       12      12

## Save the data for later use
# dir.create("rda", showWarnings = FALSE)
# save(twas, file = "rda/twas.Rdata")

## Read in the RSE info
## adapted from https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/load_funs.R
load_rse <- function(type) {

    # expmnt data
    load_file <- Sys.glob(here("dev_twas", "filter_data", paste0(opt$region, "_rda"), paste0(opt$region, "_", opt$feature, "_", "hg38_rseGene_n*.RData")))

    # load_file <- file.path(
    #     "/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/",
    #     paste0("rse_", type, "_unfiltered.Rdata")
    # )

    stopifnot(file.exists(load_file))
    load(load_file)

    ## Get the appropriate object
    rse <- rse_gene

    return(rse)
}

rse <- map(unique(twas$all$feature), load_rse)
names(rse) <- unique(twas$all$feature)

## Add the gene id
twas_exp <- map(twas, function(tw) {
    ## For testing the inner part of the function
    # tw <- twas$included

    by_feat <- split(tw, tw$feature)
    ## Make sure it's in the right order
    if (!identical(names(by_feat), names(rse))) {
        message(paste(Sys.time(), "fixing order"))
        by_feat <- by_feat[names(rse)]
    }

    ## Now add the gene gencode ID and symbol
    result <- pmap_dfr(
        list(by_feat, rse, names(rse)),
        function(info, rs, feature) {

            ## Find the appropriate variables
            gene_var <- case_when(
                feature == "gene" ~ "gencodeID",
                feature == "exon" ~ "gencodeID",
                feature == "jxn" ~ "newGeneID",
                feature == "tx" ~ "gene_id"
            )
            symbol_var <- case_when(
                feature == "gene" ~ "Symbol",
                feature == "exon" ~ "Symbol",
                feature == "jxn" ~ "newGeneSymbol",
                feature == "tx" ~ "gene_name"
            )

            ## Match by id
            m <- match(info$ID, names(rowRanges(rs)))
            stopifnot(!is.na(m))

            ## Add the gene id/symbol
            info$geneid <- mcols(rowRanges(rs))[[gene_var]][m]
            info$genesymbol <- mcols(rowRanges(rs))[[symbol_var]][m]

            ## Done
            return(info)
        }
    )

    return(result)
})
names(twas_exp) <- names(twas)

## Save the data for later use
dir.create("rda", showWarnings = FALSE)
save(twas_exp, file = "rda/twas_exp.Rdata")

################################################ end of issue #1


source("/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/twas_functions.R")

load("rda/twas_exp.Rdata", verbose = TRUE)

## Andrew's exploration code that focuses on the 'all' part
tt <- twas_exp$all
## Drop TWAS NA p-values
tt <- tt[!is.na(tt$TWAS.P), ]
## Focus on CLOZUK+PGC2 (psycm) GWAS
tt <- tt[which(tt$type == "psycm"), ]

## Add GWAS p-value and OR from the original sumstats file
original <- read_tsv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.sumstats.txt")
original$CHR[original$CHR == 23] <- "X"
original$hg19_pos <- with(original, paste0(CHR, ":", BP))

snpmap <- read_tsv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname",
    col_types = cols(
        SNP = col_character(),
        A1 = col_character(),
        A2 = col_character(),
        Z = col_double(),
        N = col_double(),
        chr = col_character(),
        basepair = col_double(),
        basepairhg19 = col_double(),
        originalSNP = col_character()
    )
)

## Load big snpMap table
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
snpMap$pos_hg38_info <- paste0(gsub("^chr", "", snpMap$chr_hg38), ":", snpMap$pos_hg38)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd <- which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap <- snpMap[-snpInd, ]

## Match to the big snpMap table
m_to_fullMap <- match(tt$BEST.GWAS.ID, snpMap$SNP)
stopifnot(!any(is.na(m_to_fullMap)))
m_to_fullMap_qtl <- match(tt$EQTL.ID, snpMap$SNP)
stopifnot(!any(is.na(m_to_fullMap_qtl)))

## Add pos hg19 and pos hg38
tt$BEST.GWAS.pos_hg19 <- snpMap$pos_hg19[m_to_fullMap]
tt$BEST.GWAS.pos_hg38 <- snpMap$pos_hg38_info[m_to_fullMap]

tt$EQTL.pos_hg19 <- snpMap$pos_hg19[m_to_fullMap_qtl]
tt$EQTL.pos_hg38 <- snpMap$pos_hg38_info[m_to_fullMap_qtl]




m_to_map <- match(tt$BEST.GWAS.ID, snpmap$SNP)
table(is.na(tt$BEST.GWAS.ID))
#  FALSE
# 723436
table(is.na(m_to_map))
## Hm... I'm not sure why some are NAs
#  FALSE   TRUE
# 717610   5826
print(tt[head(which(is.na(m_to_map))), ], width = 200)

## Hm....
m_to_map_qtl <- match(tt$EQTL.ID, snpmap$SNP)
table(is.na(m_to_map_qtl))
#  FALSE   TRUE
# 693795   7170

addmargins(table(
    "By BEST.GWAS.ID" = is.na(m_to_map),
    "By EQTL.ID" = is.na(m_to_map_qtl)
))
#                By EQTL.ID
# By BEST.GWAS.ID  FALSE   TRUE    Sum
#           FALSE 710348   7262 717610
#           TRUE    5623    203   5826
#           Sum   715971   7465 723436

## Well, after that it all looks ok
m_to_ori <- match(snpmap$originalSNP[m_to_map], original$SNP)
table(is.na(m_to_ori))
#  FALSE   TRUE
# 717610   5826

## Hm... it's odd that the same number don't match by either name or chr position
m_to_ori2 <- match(tt$BEST.GWAS.pos_hg19, original$hg19_pos)
table(is.na(m_to_ori), is.na(m_to_ori2))
#        FALSE   TRUE
# FALSE 717610      0
# TRUE       0   5826
## Supplement m_to_ori (by name) with the matching by chr and position in hg19
m_to_ori[is.na(m_to_ori)] <- m_to_ori2[is.na(m_to_ori)]
table(is.na(m_to_ori))
#  FALSE   TRUE
# 717610   5826

m_to_map_qtl2 <- match(tt$EQTL.pos_hg19, original$hg19_pos)
table(is.na(m_to_map_qtl), is.na(m_to_map_qtl2))
#        FALSE   TRUE
# FALSE 715971      0
# TRUE       0   7465



## Lets get the originally reported summarized p-values
BEST.GWAS.P <- original$P[m_to_ori]
## This is how the calculate the p-values displayed by FUSION-TWAS
# https://github.com/gusevlab/fusion_twas/blob/master/FUSION.post_process.R#L641
BEST.GWAS.P.computed <- 2 * (pnorm(abs(tt$BEST.GWAS.Z), lower.tail = F))

## Ok, assign them to our table
tt$BEST.GWAS.P <- original$P[m_to_ori]
tt$BEST.GWAS.OR <- original$OR[m_to_ori]
tt$BEST.GWAS.SE <- original$SE[m_to_ori]
tt$BEST.GWAS.P.computed <- 2 * (pnorm(abs(tt$BEST.GWAS.Z), lower.tail = F))
tt$EQTL.P.computed <- 2 * (pnorm(abs(tt$EQTL.GWAS.Z), lower.tail = F))

## Compute FDR/Bonf by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, "fdr")
        reg_feat$TWAS.Bonf <- p.adjust(reg_feat$TWAS.P, "bonf")
        reg_feat$BEST.GWAS.FDR <- p.adjust(reg_feat$BEST.GWAS.P, "fdr")
        reg_feat$BEST.GWAS.FDR.computed <- p.adjust(reg_feat$BEST.GWAS.P.computed, "fdr")
        reg_feat$EQTL.FDR.computed <- p.adjust(reg_feat$EQTL.P.computed, "fdr")
        return(reg_feat)
    })
    return(res[order(res$TWAS.P), ])
})

## Add raggr data
## Code based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R


## risk loci from PGC paper
indexLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors = FALSE)
indexLoci$hg19POS <- paste0(indexLoci$Chromosome, ":", indexLoci$snp_pos_hg19)

## risk loci from PGC paper + rAggr proxy markers
riskLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors = FALSE)
colnames(riskLoci) <- gsub("\\.", "_", colnames(riskLoci))
length(unique(riskLoci$SNP2_Name))
# [1] 10981
riskLoci$hg19POS1 <- paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos)
riskLoci$hg19POS2 <- paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos)
length(unique(riskLoci$hg19POS2))
# [1] 10975


addmargins(table(
    "proxy" = snpMap$pos_hg19 %in% riskLoci$hg19POS2,
    "index" = snpMap$pos_hg19 %in% indexLoci$hg19POS
))
#        index
# proxy     FALSE    TRUE     Sum
#   FALSE 7014124       0 7014124
#   TRUE     9600     135    9735
#   Sum   7023724     135 7023859


snpMap$Status <- "Other"
snpMap$Status[snpMap$pos_hg19 %in% riskLoci$hg19POS2] <- "Proxy"
snpMap$Status[snpMap$pos_hg19 %in% indexLoci$hg19POS] <- "Index"
table(snpMap$Status)
# Index   Other   Proxy
#   135 7014124    9600

## One of the index SNPs has 2 names
length(unique(riskLoci$SNP1_Name))
# [1] 180
which(table(gsub(":.*", "", unique(riskLoci$SNP1_Name))) > 1)
# rs1023497
#         5
## Turns out that it's multi-allelic
unique(riskLoci$SNP1_Name)[grep("rs1023497", unique(riskLoci$SNP1_Name))]
# [1] "rs1023497:42340508:C:G" "rs1023497:42340508:C:A"

tt <- as_tibble(cbind(
    tt,
    get_proxy_info(tt$BEST.GWAS.pos_hg19, "BEST.GWAS."),
    get_proxy_info(tt$EQTL.pos_hg19, "EQTL.")
))
## BEST.GWAS info
# status
#  Index  Other  Proxy
#  12641 679698  31097
# [1] 131
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -339607   -8799       0   -4774    5201  474396
# status
#  Index  Other  Proxy
#     60 720661   2715
# [1] 110
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -426653  -70824       0   -4646   39548  492315

## Label the 'other' with genome significant p-values as proxy
table(tt$BEST.GWAS.P.computed < 5e-8, tt$BEST.GWAS.status)
#        Index  Other  Proxy
# FALSE      0 668369   2403
# TRUE   12641  11329  28694
tt$BEST.GWAS.status[tt$BEST.GWAS.P.computed < 5e-8 & tt$BEST.GWAS.status == "Other"] <- "Proxy"
table(tt$BEST.GWAS.P.computed < 5e-8, tt$BEST.GWAS.status)
#        Index  Other  Proxy
# FALSE      0 668369   2403
# TRUE   12641      0  40023

print(tt, width = 600)
head(as.data.frame(tt))


table(table(tt$ID))
#      1      2
# 284154 219641

ids <- unique(tt$ID)
is_DLPFC <- tt$region == "DLPFC"

i_DLPFC <- which(is_DLPFC)[match(ids, tt$ID[is_DLPFC])]
i_HIPPO <- which(!is_DLPFC)[match(ids, tt$ID[!is_DLPFC])]
stopifnot(
    all(table(tt$region) - c(
        sum(tt$region[i_DLPFC] == "DLPFC", na.rm = TRUE),
        sum(tt$region[i_HIPPO] == "HIPPO", na.rm = TRUE)
    ) == 0)
)

ttReg_map <- data.frame(
    ID = ids,
    i_DLPFC = i_DLPFC,
    i_HIPPO = i_HIPPO,
    stringsAsFactors = FALSE
)

region_twas_z <- get_variable_by_region("TWAS.Z", NAs_0 = TRUE)

## How I figured out something weird :P
## Some values on the cross with low/high Zs were labeled as "None":
## they were NAs on one of the two regions
head(subset(region_twas_z, FDR.5perc == "None" & DLPFC < -5))


table(region_twas_z$in_both, useNA = "ifany")
#  FALSE   TRUE
# 284154 219641
table(region_twas_z$in_both) / nrow(region_twas_z) * 100
#   FALSE    TRUE
# 56.4027 43.5973
addmargins(table(
    "in both" = region_twas_z$in_both,
    "FDR < 0.05" = region_twas_z$FDR.5perc,
    useNA = "ifany"
))
#        FDR < 0.05
# in both   None  DLPFC  HIPPO   Both    Sum
#   FALSE 264827  10587   8740      0 284154
#   TRUE  195241   8766   8264   7370 219641
#   Sum   460068  19353  17004   7370 503795
addmargins(table(
    "in both" = region_twas_z$in_both,
    "Bonf < 0.05" = region_twas_z$Bonf.5perc,
    useNA = "ifany"
))
#        Bonf < 0.05
# in both   None  DLPFC  HIPPO   Both    Sum
#   FALSE 282343    997    814      0 284154
#   TRUE  217027    978    939    697 219641
#   Sum   499370   1975   1753    697 503795

table(region_twas_z$BEST.GWAS.status)
#  Other Risk Locus
# 465196      38599


dir.create("pdf", showWarnings = FALSE)
pdf("pdf/twas_z.pdf", useDingbats = FALSE, width = 24, height = 14)
ggplot(
    region_twas_z,
    aes(x = DLPFC, y = HIPPO, color = FDR.5perc, shape = in_both)
) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle("TWAS Z by brain region") +
    scale_color_manual(values = c("grey80", "dark orange", "skyblue3", "purple"))

ggplot(
    region_twas_z,
    aes(x = DLPFC, y = HIPPO, color = Bonf.5perc, shape = in_both)
) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle("TWAS Z by brain region") +
    scale_color_manual(values = c("grey80", "dark orange", "skyblue3", "purple"))
dev.off()

pdf("pdf/twas_z_gene.pdf", useDingbats = FALSE, width = 10, height = 10)
ggplot(
    subset(region_twas_z, feature == "gene"),
    aes(x = DLPFC, y = HIPPO, color = FDR.5perc, shape = in_both)
) +
    geom_point() +
    # facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle("TWAS Z by brain region") +
    scale_color_manual(values = c("grey80", "dark orange", "skyblue3", "purple"))

ggplot(
    subset(region_twas_z, feature == "gene"),
    aes(x = DLPFC, y = HIPPO, color = Bonf.5perc, shape = in_both)
) +
    geom_point() +
    # facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle("TWAS Z by brain region") +
    scale_color_manual(values = c("grey80", "dark orange", "skyblue3", "purple"))
dev.off()

## Find the discordant ones
options(width = 300)
subset(region_twas_z, FDR.5perc == "Both" & sign(DLPFC) != sign(HIPPO))
## Lots of output
dim(subset(region_twas_z, FDR.5perc == "Both" & sign(DLPFC) != sign(HIPPO)))
# [1] 338  16

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=LINC01378&keywords=ENSG00000236922
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=TMEM98&keywords=ENSG00000006042

## Get the numbers of points in the different parts of the plot
map(
    split(region_twas_z, region_twas_z$feature),
    ~ map(
        split(.x, .x$BEST.GWAS.status),
        ~ addmargins(table("FDR <5%" = .x$FDR.5perc, "In both" = .x$in_both))
    )
)
# $gene
# $gene$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   6379  9229 15608
#   DLPFC   168   266   434
#   HIPPO   140   284   424
#   Both      0   237   237
#   Sum    6687 10016 16703
#
# $gene$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None    384  447  831
#   DLPFC    88   87  175
#   HIPPO    66   82  148
#   Both      0  142  142
#   Sum     538  758 1296
#
#
# $exon
# $exon$Other
#        In both
# FDR <5%  FALSE   TRUE    Sum
#   None  124311  97652 221963
#   DLPFC   3664   3543   7207
#   HIPPO   2895   3401   6296
#   Both       0   2213   2213
#   Sum   130870 106809 237679
#
# $exon$`Risk Locus`
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   7989  5074 13063
#   DLPFC  1845  1191  3036
#   HIPPO  1505  1083  2588
#   Both      0  1625  1625
#   Sum   11339  8973 20312
#
#
# $jxn
# $jxn$Other
#        In both
# FDR <5%  FALSE   TRUE    Sum
#   None   90117  57411 147528
#   DLPFC   2497   2083   4580
#   HIPPO   2089   1889   3978
#   Both       0   1324   1324
#   Sum    94703  62707 157410
#
# $jxn$`Risk Locus`
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   5495  2774  8269
#   DLPFC  1188   583  1771
#   HIPPO   975   604  1579
#   Both      0   863   863
#   Sum    7658  4824 12482
#
#
# $tx
# $tx$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  28286 21557 49843
#   DLPFC   733   777  1510
#   HIPPO   744   724  1468
#   Both      0   583   583
#   Sum   29763 23641 53404
#
# $tx$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None   1866 1097 2963
#   DLPFC   404  236  640
#   HIPPO   326  197  523
#   Both      0  383  383
#   Sum    2596 1913 4509



## Now for Bonf
map(
    split(region_twas_z, region_twas_z$feature),
    ~ map(
        split(.x, .x$BEST.GWAS.status),
        ~ addmargins(table("Bonf <5%" = .x$Bonf.5perc, "In both" = .x$in_both))
    )
)
# $gene
# $gene$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None   6672  9974 16646
#    DLPFC     6    18    24
#    HIPPO     9    17    26
#    Both      0     7     7
#    Sum    6687 10016 16703
#
# $gene$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None    480  627 1107
#    DLPFC    33   37   70
#    HIPPO    25   40   65
#    Both      0   54   54
#    Sum     538  758 1296
#
#
# $exon
# $exon$Other
#         In both
# Bonf <5%  FALSE   TRUE    Sum
#    None  130784 106687 237471
#    DLPFC     48     56    104
#    HIPPO     38     62    100
#    Both       0      4      4
#    Sum   130870 106809 237679
#
# $exon$`Risk Locus`
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  10570  7747 18317
#    DLPFC   429   459   888
#    HIPPO   340   440   780
#    Both      0   327   327
#    Sum   11339  8973 20312
#
#
# $jxn
# $jxn$Other
#         In both
# Bonf <5%  FALSE   TRUE    Sum
#    None   94639  62627 157266
#    DLPFC     35     42     77
#    HIPPO     29     33     62
#    Both       0      5      5
#    Sum    94703  62707 157410
#
# $jxn$`Risk Locus`
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None   7099  4168 11267
#    DLPFC   300   249   549
#    HIPPO   259   224   483
#    Both      0   183   183
#    Sum    7658  4824 12482
#
#
# $tx
# $tx$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  29721 23589 53310
#    DLPFC    22    19    41
#    HIPPO    20    25    45
#    Both      0     8     8
#    Sum   29763 23641 53404
#
# $tx$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None   2378 1608 3986
#    DLPFC   124   98  222
#    HIPPO    94   98  192
#    Both      0  109  109
#    Sum    2596 1913 4509

corrs <- map(split(region_twas_z, factor(region_twas_z$feature, levels = features)), ~ map_dbl(split(.x, ifelse(.x$BEST.GWAS.status_DLPFC == "Other", "Other", "Risk Loci")), ~ cor(.x$DLPFC[.x$in_both], .x$HIPPO[.x$in_both])))
corrs
# $gene
#     Other Risk Loci
# 0.6379117 0.6613453
#
# $exon
#     Other Risk Loci
# 0.5398872 0.5834149
#
# $jxn
#     Other Risk Loci
# 0.5194039 0.5787448
#
# $tx
#     Other Risk Loci
# 0.5668870 0.6427932

range(unlist(corrs))
# [1] 0.5194039 0.6613453

## Load SCZD case-control data
## From https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_case_control.R#L65-L71
outFeat <- lapply(c("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda", "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda"), function(f) {
    message(paste(Sys.time(), "loading", f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub("\\..*", "", outTx$gene_id)
    return(list("gene" = outGene, "exon" = outExon, "jxn" = outJxn, "tx" = outTx))
})
names(outFeat) <- c("DLPFC", "HIPPO")

## Moved this up so it'll be saved in the object
tt$status <- ifelse(tt$BEST.GWAS.status == "Other", "Other", "Risk Locus")


## Add the SCZD info back to the tt object
tt$SCZD_FDR <- tt$SCZD_pvalue <- tt$SCZD_t <- NA
for (feat in features) {
    for (region in names(outFeat)) {
        message(paste(Sys.time(), "processing", region, "at the", feat, "level"))
        i <- which(tt$feature == feat & tt$region == region)
        current <- outFeat[[region]][[feat]]
        print(length(i))
        # table(tt$region[i], tt$feature[i])
        m <- match(tt$ID[i], rownames(current))
        stopifnot(!any(is.na(m)))
        tt$SCZD_t[i] <- current$t[m]
        tt$SCZD_pvalue[i] <- current$P.Value[m]
        tt$SCZD_FDR[i] <- current$adj.P.Val[m]
    }
}
## For checking the code
# sapply(tt[, 49:51], function(x) { sum(is.na(x)) })


## Subset by significant (TWAS FDR < 5%)
ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map_dfr(ttSig, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1 26723 24374
# 2    52    52

ttSig_bonf <- map(split(tt, tt$region), ~ .x[.x$TWAS.Bonf < 0.05, ])
map_dfr(ttSig_bonf, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1  2672  2450
# 2    52    52

map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  5837  5754

map_int(ttSig_bonf, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#   682   696

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   433   465

dim(tt)
# [1] 723436     52
dim(region_twas_z)
# [1] 503795     16


## save for later
dir.create("rda", showWarnings = FALSE)
save(tt, ttSig, ttSig_bonf, get_variable_by_region, ttReg_map, region_twas_z, file = "rda/tt_objects.Rdata")






tt_sigonly <- tt[tt$TWAS.FDR < 0.05, ]
with(tt_sigonly, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = "ifany")))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy   Sum
#            Index    48  3615   563  4226
#            Other     0 34598    13 34611
#            Proxy     7 10332  1921 12260
#            Sum      55 48545  2497 51097

tt_sigonly_bonf <- tt[tt$TWAS.Bonf < 0.05, ]
with(tt_sigonly_bonf, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = "ifany")))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index    34   838   288 1160
#            Other     0   527     0  527
#            Proxy     7  2152  1276 3435
#            Sum      41  3517  1564 5122

create_gwas_or_eqtl(tt_sigonly, "pdf/twas_fdr5perc_vs_gwas_or_eqtl.pdf", "FDR")
create_gwas_or_eqtl(tt_sigonly_bonf, "pdf/twas_bonf5perc_vs_gwas_or_eqtl.pdf", "Bonf")

create_by_status(tt_sigonly, "pdf/twas_fdr5perc_by_status.pdf", "FDR")
create_by_status(tt_sigonly_bonf, "pdf/twas_bonf5perc_by_status.pdf", "Bonf")


map(ttSig, ~ map(split(.x, .x$feature), ~
addmargins(table(
    "BEST GWAS P < 5e-08" = .x$BEST.GWAS.P.computed < 5e-08,
    "Risk Locus (by BEST GWAS)" =  .x$BEST.GWAS.status != "Other",
    useNA = "ifany"
))))
# $DLPFC
# $DLPFC$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE  TRUE   Sum
#               FALSE  9420   127  9547
#               TRUE      0  4534  4534
#               Sum    9420  4661 14081
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   671   13 684
#               TRUE      0  304 304
#               Sum     671  317 988
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  5904   86 5990
#               TRUE      0 2548 2548
#               Sum    5904 2634 8538
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  2093   25 2118
#               TRUE      0  998  998
#               Sum    2093 1023 3116
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE  TRUE   Sum
#               FALSE  8509   108  8617
#               TRUE      0  4105  4105
#               Sum    8509  4213 12722
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   661   12 673
#               TRUE      0  278 278
#               Sum     661  290 951
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  5302   82 5384
#               TRUE      0 2360 2360
#               Sum    5302 2442 7744
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  2051   35 2086
#               TRUE      0  871  871
#               Sum    2051  906 2957


map(ttSig_bonf, ~ map(split(.x, .x$feature), ~
addmargins(table(
    "BEST GWAS P < 5e-08" = .x$BEST.GWAS.P.computed < 5e-08,
    "Risk Locus (by BEST GWAS)" =  .x$BEST.GWAS.status != "Other",
    useNA = "ifany"
))))
# $DLPFC
# $DLPFC$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   108    2  110
#               TRUE      0 1213 1213
#               Sum     108 1215 1323
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    31    1  32
#               TRUE      0  123 123
#               Sum      31  124 155
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    82    5  87
#               TRUE      0  727 727
#               Sum      82  732 814
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    49    3  52
#               TRUE      0  328 328
#               Sum      49  331 380
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   104    1  105
#               TRUE      0 1106 1106
#               Sum     104 1107 1211
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    33    2  35
#               TRUE      0  117 117
#               Sum      33  119 152
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    67    4  71
#               TRUE      0  662 662
#               Sum      67  666 733
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    53    2  55
#               TRUE      0  299 299
#               Sum      53  301 354


## Add locus considered section

## Read in the files that Emily cleaned up at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R
raggr_clean_files <- c(
    "HIPPO" = "/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_hippo_eqtls_fdr01.csv",
    "DLPFC" = "/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_dlp_eqtls_fdr01.csv"
)
raggr_clean <- map(raggr_clean_files, read.csv, stringsAsFactors = FALSE)
names(raggr_clean) <- names(raggr_clean_files)


clean_by_state(
    map2_dfr(
        raggr_clean,
        split(tt, factor(tt$region, levels = names(raggr_clean))),
        ~ table(unique(.x$IndexSNP) %in% .y$BEST.GWAS.indexSNP)
    )
)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    17    23 FALSE
# 2    86    93 TRUE

stopifnot(identical(by_locus(1.1), by_locus(1.1, var = "TWAS.Bonf")))


## raggr eQTL locus considered in the TWAS analysis
by_locus(1.1)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    17    23 FALSE
# 2    86    93 TRUE
perc_locus(1.1)
#      HIPPO    DLPFC
# 1 83.49515 80.17241

## raggr eQTL locus that have a TWAS FDR <5% result
by_locus(0.05)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    24    28 FALSE
# 2    79    88 TRUE
perc_locus(0.05)
#      HIPPO    DLPFC
# 1 76.69903 75.86207
perc_locus(0.05) / perc_locus(1.1) * 100
#      HIPPO    DLPFC
# 1 91.86047 93.54839

## raggr eQTL locus that have a TWAS Bonf <5% result
by_locus(0.05, var = "TWAS.Bonf")
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    42    48 FALSE
# 2    61    68 TRUE
perc_locus(0.05, var = "TWAS.Bonf")
#     HIPPO    DLPFC
# 1 59.2233 58.62069
perc_locus(0.05, var = "TWAS.Bonf") / perc_locus(1.1) * 100
#      HIPPO    DLPFC
# 1 70.93023 73.11828


## Overlap by region of all loci considered
shared_by_locus(1.1)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   8     1     0
# 11  78     1     1

## Now with the ones that are TWAS FDR < 5%
shared_by_locus(0.05)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   6     1     0
# 11  73     1     1

## Now with the ones that are TWAS Bonf < 5%
shared_by_locus(0.05, var = "TWAS.Bonf")
#    num HIPPO DLPFC
# 00   0     0     0
# 01  14     0     1
# 10   7     1     0
# 11  54     1     1


pdf("pdf/venn_by_locus.pdf", useDingbats = FALSE)
make_pretty_venn(1.1, "rAggr loci considered in TWAS")
make_pretty_venn(0.05, "rAggr loci with TWAS FDR<5%")
make_pretty_venn(0.05, "rAggr loci with TWAS Bonf<5%", var = "TWAS.Bonf")
dev.off()
system("rm VennDiagram*")

## Features by locus that were considered in the TWAS analysis
clean_tabs(by_feature(cut = 1.1))
## Numbers differ from https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas.R#L774-L784
## because in this script we already dropped the results with NA TWAS.P values
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1     2     2 FALSE gene
# 2    48    65 TRUE  gene
# 3     2     1 FALSE exon
# 4    75    86 TRUE  exon
# 5     7     5 FALSE jxn
# 6    81    93 TRUE  jxn
# 7     3     1 FALSE tx
# 8    62    75 TRUE  tx

perc_feature(1.1)
#      HIPPO    DLPFC feature
# 1 96.00000 97.01493    gene
# 2 97.40260 98.85057    exon
# 3 92.04545 94.89796     jxn
# 4 95.38462 98.68421      tx

## Features by locus that have a TWAS FDR<5% result
clean_tabs(by_feature(cut = 0.05))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    14    14 FALSE gene
# 2    36    53 TRUE  gene
# 3    13    13 FALSE exon
# 4    64    74 TRUE  exon
# 5    15    15 FALSE jxn
# 6    73    83 TRUE  jxn
# 7    15    11 FALSE tx
# 8    50    65 TRUE  tx

perc_feature(0.05)
#      HIPPO    DLPFC feature
# 1 72.00000 79.10448    gene
# 2 83.11688 85.05747    exon
# 3 82.95455 84.69388     jxn
# 4 76.92308 85.52632      tx

## Compute the percent using as denominator the number of features considered
cbind(perc_feature(0.05)[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 75.00000 81.53846    gene
# 2 85.33333 86.04651    exon
# 3 90.12346 89.24731     jxn
# 4 80.64516 86.66667      tx

## Now with Bonf < 5%
clean_tabs(by_feature(cut = 0.05, var = "TWAS.Bonf"))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    25    31 FALSE gene
# 2    25    36 TRUE  gene
# 3    38    34 FALSE exon
# 4    39    53 TRUE  exon
# 5    41    41 FALSE jxn
# 6    47    57 TRUE  jxn
# 7    28    30 FALSE tx
# 8    37    46 TRUE  tx

perc_feature(0.05, var = "TWAS.Bonf")
#      HIPPO    DLPFC feature
# 1 50.00000 53.73134    gene
# 2 50.64935 60.91954    exon
# 3 53.40909 58.16327     jxn
# 4 56.92308 60.52632      tx

cbind(perc_feature(0.05, var = "TWAS.Bonf")[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 52.08333 55.38462    gene
# 2 52.00000 61.62791    exon
# 3 58.02469 61.29032     jxn
# 4 59.67742 61.33333      tx

## Find the locus that have shared features
## between rAggr and TWAS (given a FDR cut)
## then find the locus overlap across regions.
shared_by_locus_by_feature(1.1)
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  22     0     1
# 10   5     1     0
# 11  43     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  20     0     1
# 10   9     1     0
# 11  66     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  22     0     1
# 10  10     1     0
# 11  71     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  23     0     1
# 10  10     1     0
# 11  52     1     1

## TWAS FDR <5%
shared_by_locus_by_feature(0.05)
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  19     0     1
# 10   2     1     0
# 11  34     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  19     0     1
# 10   9     1     0
# 11  55     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  18     0     1
# 10   8     1     0
# 11  65     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  21     0     1
# 10   6     1     0
# 11  44     1     1

## TWAS Bonf <5%
shared_by_locus_by_feature(0.05, var = "TWAS.Bonf")
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  14     0     1
# 10   3     1     0
# 11  22     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  18     0     1
# 10   4     1     0
# 11  35     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  13     0     1
# 10   3     1     0
# 11  44     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  18     0     1
# 10   9     1     0
# 11  28     1     1


pdf("pdf/venn_by_locus_by_feature.pdf", useDingbats = FALSE)
make_pretty_venn_by_feature(1.1, "rAggr loci considered in TWAS")
make_pretty_venn_by_feature(0.05, "rAggr loci with TWAS FDR<5%")
make_pretty_venn_by_feature(0.05, "rAggr loci with TWAS Bonf<5%", var = "TWAS.Bonf")
dev.off()
system("rm VennDiagram*")


## Venn diagrams of features by region, then joint (grouped by gene id)

## Number of features that have TWAS weights
n_feat_considered <- cbind(map_dfr(split(tt, tt$region), ~
map_dbl(
    split(.x, factor(.x$feature, levels = features)),
    ~ nrow(.x)
)), features)
n_feat_considered
#    DLPFC  HIPPO features
# 1  14946  13827     gene
# 2 193190 180583     exon
# 3 121949 115474      jxn
# 4  41847  41620       tx

## Percent from total number of features considered
# $ grep -A 1 "Final RSE" logs/compute_weights_full_DLPFC_*
# compute_weights_DLPFC_exon.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_exon.txt-[1] 381196    397
# --
# compute_weights_DLPFC_gene.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_gene.txt-[1] 23402   397
# --
# compute_weights_DLPFC_jxn.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_jxn.txt-[1] 269261    397
# --
# compute_weights_DLPFC_tx.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_tx.txt-[1] 88969   397

## In percent of features passed to the weight computation step
n_feat_total <- c("gene" = 23402, "exon" = 381196, "jxn" = 269261, "tx" = 88969)
cbind(map_dfr(n_feat_considered[, 1:2], ~ .x / n_feat_total * 100), features)
#      DLPFC    HIPPO features
# 1 63.86634 59.08469     gene
# 2 50.67997 47.37274     exon
# 3 45.29026 42.88553      jxn
# 4 47.03548 46.78034       tx

## Number of features with TWAS FDR<5%
n_feat_twas5perc <- cbind(map_dfr(ttSig, ~
map_dbl(
    split(.x, factor(.x$feature, levels = features)),
    ~ nrow(.x)
)), features)
n_feat_twas5perc
#   DLPFC HIPPO features
# 1   988   951     gene
# 2 14081 12722     exon
# 3  8538  7744      jxn
# 4  3116  2957       tx

cbind(n_feat_twas5perc[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#      DLPFC    HIPPO features
# 1 6.610464 6.877848     gene
# 2 7.288680 7.044960     exon
# 3 7.001287 6.706272      jxn
# 4 7.446173 7.104757       tx

## Now with TWAS Bonf <5%
n_feat_twas5perc_bonf <- cbind(map_dfr(ttSig_bonf, ~
map_dbl(
    split(.x, factor(.x$feature, levels = features)),
    ~ nrow(.x)
)), features)
n_feat_twas5perc_bonf
#   DLPFC HIPPO features
# 1   155   152     gene
# 2  1323  1211     exon
# 3   814   733      jxn
# 4   380   354       tx

cbind(n_feat_twas5perc_bonf[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#       DLPFC     HIPPO features
# 1 1.0370668 1.0992985     gene
# 2 0.6848181 0.6706058     exon
# 3 0.6674921 0.6347749      jxn
# 4 0.9080699 0.8505526       tx

shared_by_feature(1.1)
# $gene
#      num DLPFC HIPPO
# 00     0     0     0
# 01  3053     0     1
# 10  4172     1     0
# 11 10774     1     1
#
# $exon
#       num DLPFC HIPPO
# 00      0     0     0
# 01  64801     0     1
# 10  77408     1     0
# 11 115782     1     1
#
# $jxn
#      num DLPFC HIPPO
# 00     0     0     0
# 01 47943     0     1
# 10 54418     1     0
# 11 67531     1     1
#
# $tx
#      num DLPFC HIPPO
# 00     0     0     0
# 01 16066     0     1
# 10 16293     1     0
# 11 25554     1     1

## TWAS FDR <5%
shared_by_feature(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01 572     0     1
# 10 609     1     0
# 11 379     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01  8884     0     1
# 10 10243     1     0
# 11  3838     1     1
#
# $jxn
#     num DLPFC HIPPO
# 00    0     0     0
# 01 5557     0     1
# 10 6351     1     0
# 11 2187     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1991     0     1
# 10 2150     1     0
# 11  966     1     1

## TWAS Bonf <5%
shared_by_feature(0.05, var = "TWAS.Bonf")
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  91     0     1
# 10  94     1     0
# 11  61     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01 880     0     1
# 10 992     1     0
# 11 331     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 545     0     1
# 10 626     1     0
# 11 188     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 237     0     1
# 10 263     1     0
# 11 117     1     1



pdf("pdf/venn_by_feature.pdf", useDingbats = FALSE)
make_pretty_venn_shared_by_feature(1.1, "Features with TWAS weights")
make_pretty_venn_shared_by_feature(0.05, "Features with TWAS FDR<5%")
make_pretty_venn_shared_by_feature(0.05, "Features with TWAS Bonf<5%", var = "TWAS.Bonf")
dev.off()
system("rm VennDiagram*")


## Now by gene id
shared_by_geneid(1.1)
# $gene
#      num DLPFC HIPPO
# 00     0     0     0
# 01  3053     0     1
# 10  4172     1     0
# 11 10774     1     1
#
# $exon
#       num DLPFC HIPPO
# 00      0     0     0
# 01   2570     0     1
# 10   3001     1     0
# 11 190189     1     1
#
# $jxn
#       num DLPFC HIPPO
# 00      0     0     0
# 01   1758     0     1
# 10   1769     1     0
# 11 107019     1     1
#
# $tx
#      num DLPFC HIPPO
# 00     0     0     0
# 01  4652     0     1
# 10  4462     1     0
# 11 37385     1     1

## TWAS FDR <5%
shared_by_geneid(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01 572     0     1
# 10 609     1     0
# 11 379     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01  2799     0     1
# 10  3297     1     0
# 11 10784     1     1
#
# $jxn
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1669     0     1
# 10 1832     1     0
# 11 5866     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1276     0     1
# 10 1349     1     0
# 11 1767     1     1

## TWAS Bonf <5%
shared_by_geneid(0.05, var = "TWAS.Bonf")
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  91     0     1
# 10  94     1     0
# 11  61     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01 284     0     1
# 10 346     1     0
# 11 977     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 192     0     1
# 10 225     1     0
# 11 523     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 156     0     1
# 10 165     1     0
# 11 215     1     1

pdf("pdf/venn_by_feature_using_geneid.pdf", useDingbats = FALSE)
make_pretty_venn_shared_by_geneid(1.1, "Features with TWAS weights (by gene ID)")
make_pretty_venn_shared_by_geneid(0.05, "Features with TWAS FDR<5% (by gene ID)")
make_pretty_venn_shared_by_geneid(0.05, "Features with TWAS Bonf<5% (by gene ID)", var = "TWAS.Bonf")
dev.off()
system("rm VennDiagram*")


pdf("pdf/venn_by_feature_using_geneid_across_features.pdf", useDingbats = FALSE)
make_pretty_venn_shared_by_geneid2(1.1, "Features with TWAS weights (by gene ID)")
make_pretty_venn_shared_by_geneid2(0.05, "Features with TWAS FDR<5% (by gene ID)")
make_pretty_venn_shared_by_geneid2(0.05, "Features with TWAS Bonf<5% (by gene ID)", var = "TWAS.Bonf")
dev.off()
system("rm VennDiagram*")

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.FDR < 0.05))), region = c("DLPFC", "HIPPO"))
#    exon gene  jxn   tx region
# 1 14081  988 8538 3116  DLPFC
# 2 12722  951 7744 2957  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.P < 5e-08))), region = c("DLPFC", "HIPPO"))
#   exon gene jxn  tx region
# 1  974   64 553 213  DLPFC
# 2  876   68 498 190  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(p.adjust(.x$TWAS.P, "bonf") < 0.05))), region = c("DLPFC", "HIPPO"))
#   exon gene jxn  tx region
# 1 1323  155 814 380  DLPFC
# 2 1211  152 733 354  HIPPO



## Compare SCZD t vs TWAS z
pdf("pdf/sczd_t_vs_twas_z.pdf", useDingbats = FALSE, width = 21, height = 18)

## Used https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
dat_text <- map_dfr(features, function(feat) {
    map_dfr(names(outFeat), function(region) {
        map_dfr(c("Other", "Risk Locus"), function(state) {
            i <- which(tt$feature == feat & tt$region == region & tt$status == state)
            data.frame(
                region = region,
                feature = feat,
                status = state,
                cor = signif(cor(tt$TWAS.Z[i], tt$SCZD_t[i]), 3),
                stringsAsFactors = FALSE
            )
        })
    })
})

## Used https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals/39625148
ggplot(tt, aes(x = TWAS.Z, y = SCZD_t, color = BEST.GWAS.P.computed < 5e-08)) +
    geom_point() +
    facet_grid(region *
        status ~
    factor(feature, levels = c("gene", "exon", "jxn", "tx"))) +
    theme_bw(base_size = 30) +
    ggtitle("TWAS vs SCZD differential expression") +
    xlab("TWAS Z score") +
    ylab("SCZD vs control t-statistic") +
    labs(caption = "Risk Loci by BEST GWAS") +
    guides(color = guide_legend(title = "BEST GWAS\np < 5e-08")) +
    geom_text(
        data = dat_text,
        color = "black",
        size = 7,
        mapping = aes(x = -4, y = 5.5, label = paste0("rho=", formatC(cor, format = "e", digits = 2))),
    )
dev.off()

pdf("pdf/sczd_t_vs_twas_z_gene.pdf", useDingbats = FALSE, width = 16, height = 14)
## Used https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals/39625148
ggplot(subset(tt, feature == "gene"), aes(x = TWAS.Z, y = SCZD_t, color = BEST.GWAS.P.computed < 5e-08)) +
    geom_point() +
    # facet_grid(region *
    #     status ~
    #     factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
    # ) +
    facet_grid(
        status ~
        region
    ) +
    theme_bw(base_size = 30) +
    ggtitle("TWAS vs SCZD differential expression") +
    xlab("TWAS Z score") +
    ylab("SCZD vs control t-statistic") +
    labs(caption = "Risk Loci by BEST GWAS") +
    guides(color = guide_legend(title = "BEST GWAS\np < 5e-08")) +
    geom_text(
        data = subset(dat_text, feature == "gene"),
        color = "black",
        size = 10,
        mapping = aes(x = -4, y = 5.5, label = paste0("rho=", formatC(cor, format = "e", digits = 2))),
    )
dev.off()

## Create the summary twas table
## Complements venn_by_feature*.pdf files
options(width = 300)
(summary_twas_tab <- make_summary_twas_table(tt))
#     type          set features genes twas_fdr_features twas_fdr_in_features twas_fdr_out_features twas_fdr_genes twas_fdr_in_genes twas_fdr_out_genes twas_bonf_features twas_bonf_in_features twas_bonf_out_features twas_bonf_genes twas_bonf_in_genes twas_bonf_out_genes
# 1   gene        DLPFC    14946 14946               988                  317                   671            988               317                671                155                   124                     31             155                124                  31
# 2   exon        DLPFC   193190 18327             14081                 4661                  9420           3751               838               2932               1323                  1215                    108             370                319                  51
# 3    jxn        DLPFC   121949 15660              8538                 2634                  5904           3236               706               2540                814                   732                     82             360                295                  65
# 4     tx        DLPFC    41847 18573              3116                 1023                  2093           2171               606               1567                380                   331                     49             269                223                  46
# 5  total        DLPFC   371932 25096             26723                 8635                 18088           5836              1151               4715               2672                  2402                    270             681                532                 151
# 6   gene        HIPPO    13827 13827               951                  290                   661            951               290                661                152                   119                     33             152                119                  33
# 7   exon        HIPPO   180583 18151             12722                 4213                  8509           3682               823               2878               1211                  1107                    104             384                326                  59
# 8    jxn        HIPPO   115474 15629              7744                 2442                  5302           3164               722               2453                733                   666                     67             343                291                  52
# 9     tx        HIPPO    41620 18727              2957                  906                  2051           2113               564               1552                354                   301                     53             260                212                  48
# 10 total        HIPPO   351504 25226             24374                 7851                 16523           5753              1153               4631               2450                  2193                    257             695                541                 155
# 11  gene Intersection    10774 10774               379                  142                   237            379               142                237                 61                    54                      7              61                 54                   7
# 12  exon Intersection   115782 16512              3838                 1625                  2213           2151               641               1519                331                   327                      4             209                191                  18
# 13   jxn Intersection    67531 14467              2187                  863                  1324           1848               538               1312                188                   183                      5             188                176                  12
# 14    tx Intersection    25554 14610               966                  383                   583            958               352                607                117                   109                      8             120                108                  12
# 15 total Intersection   219641 21144              7370                 3013                  4357           3053               819               2248                697                   673                     24             337                298                  39
# 16  gene        Union    17999 17999              1560                  465                  1095           1560               465               1095                246                   189                     57             246                189                  57
# 17  exon        Union   257991 19966             22965                 7249                 15716           5282              1020               4291               2203                  1995                    208             545                454                  92
# 18   jxn        Union   169892 16822             14095                 4213                  9882           4552               890               3681               1359                  1215                    144             515                410                 105
# 19    tx        Union    57913 22690              5107                 1546                  3561           3326               818               2512                617                   523                     94             409                327                  82
# 20 total        Union   503795 28630             43727                13473                 30254           7941              1397               6592               4425                  3922                    503             970                715                 257

## Check numbers
rbind(
    "DLPFC + HIPPO - intersection" = summary_twas_tab[5, -(1:2)] + summary_twas_tab[10, -(1:2)] - summary_twas_tab[15, -(1:2)],
    "Union" = summary_twas_tab[20, -(1:2)]
)
## Genes don't add up because those can show up more than once
## for the exon/jxn/tx features
#                              features genes twas_fdr_features twas_fdr_in_features twas_fdr_out_features twas_fdr_genes twas_fdr_in_genes twas_fdr_out_genes twas_bonf_features twas_bonf_in_features twas_bonf_out_features twas_bonf_genes twas_bonf_in_genes twas_bonf_out_genes
# DLPFC + HIPPO - intersection   503795 29178             43727                13473                 30254           8536              1485               7098               4425                  3922                    503            1039                775                 267
# Union                          503795 28630             43727                13473                 30254           7941              1397               6592               4425                  3922                    503             970                715                 257


## Save results
dir.create("summary_results", showWarnings = FALSE)
write.table(summary_twas_tab, file = "summary_results/summary_twas_table.txt", row.names = FALSE, quote = FALSE, sep = "\t")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-09-06 r77160)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-10-07
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.44.0    2019-05-02 [2] Bioconductor
#  BiocGenerics         * 0.30.0    2019-05-02 [1] Bioconductor
#  BiocParallel         * 1.18.1    2019-08-06 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  caTools                1.17.1.2  2019-03-06 [2] CRAN (R 3.6.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-09-26 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.10.0    2019-05-02 [2] Bioconductor
#  digest                 0.6.21    2019-09-20 [1] CRAN (R 3.6.1)
#  dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  formatR                1.7       2019-06-11 [1] CRAN (R 3.6.1)
#  futile.logger        * 1.4.3     2016-07-10 [1] CRAN (R 3.6.1)
#  futile.options         1.0.1     2018-04-20 [2] CRAN (R 3.6.1)
#  gdata                  2.18.0    2017-06-06 [2] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.20.0    2019-05-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.1     2019-09-09 [2] Bioconductor
#  GenomicRanges        * 1.36.1    2019-09-06 [1] Bioconductor
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  gplots               * 3.0.1.1   2019-01-27 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  gtools                 3.8.1     2018-06-26 [2] CRAN (R 3.6.1)
#  hms                    0.5.1     2019-08-23 [2] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5       2019-10-04 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.18.3    2019-09-24 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  KernSmooth             2.23-15   2015-06-29 [3] CRAN (R 3.6.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
#  lambda.r               1.2.4     2019-09-18 [1] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                * 0.3.2     2019-03-15 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  RColorBrewer         * 1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.2     2019-07-25 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readr                * 1.3.1     2018-12-21 [1] CRAN (R 3.6.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.6.1)
#  rlang                  0.4.0     2019-06-25 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-09-26 [1] Github (cloudyr/rmote@fbce611)
#  S4Vectors            * 0.22.1    2019-09-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.6.1)
#  stringr              * 1.4.0     2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.14.1    2019-07-31 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
#  vctrs                  0.2.0     2019-07-05 [1] CRAN (R 3.6.1)
#  VennDiagram          * 1.6.20    2018-03-28 [1] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.10      2019-10-01 [1] CRAN (R 3.6.1)
#  XVector                0.24.0    2019-05-02 [1] Bioconductor
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.30.0    2019-05-02 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/library
