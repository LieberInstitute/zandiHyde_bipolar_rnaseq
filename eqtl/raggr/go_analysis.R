library('SummarizedExperiment')
library('data.table')
library('clusterProfiler')
library('org.Hs.eg.db')
library('GenomicRanges')
library('sessioninfo')

dir.create('go_files', showWarnings = FALSE)

## Load the rse files
rse_files <- dir('/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/',
    pattern = 'eQTL_expressed', full.names = TRUE)
names(rse_files) <- gsub(
    'eQTL_expressed_rse_|.rda',
    '',
    dir('/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/',
        pattern = 'eQTL_expressed')
)

## rse_files have already been standarized, so I only need to load 1
load(rse_files['amygdala'], verbose = TRUE)
uni_expressed_list <- list(
    gene = unique(rowRanges(rse_gene)$ensemblID),
    exon = unique(rowRanges(rse_exon)$ensemblID),
    jxn = unique(gsub('\\..*', '', rowRanges(rse_jxn)$newGeneID)),
    tx = unique(gsub('\\..*', '', rowRanges(rse_tx)$gene_id))
)
sapply(uni_expressed_list, length)
 # gene  exon   jxn    tx
# 25136 23628 19185 27059
uni_expressed <- unique(unlist(uni_expressed_list))
length(uni_expressed)
# [1] 33932

## Just checking with DLPFC
load(rse_files['dlpfc'], verbose = TRUE)
uni_expressed_list_dlpfc <- list(
    gene = unique(rowRanges(rse_gene)$ensemblID),
    exon = unique(rowRanges(rse_exon)$ensemblID),
    jxn = unique(gsub('\\..*', '', rowRanges(rse_jxn)$newGeneID)),
    tx = unique(gsub('\\..*', '', rowRanges(rse_tx)$gene_id))
)
sapply(uni_expressed_list_dlpfc, length)
#  gene  exon   jxn    tx
# 25136 23624 19113 27059
identical(sapply(uni_expressed_list, length), sapply(uni_expressed_list_dlpfc, length))
# [1] FALSE
uni_expressed_dlpfc <- unique(unlist(uni_expressed_list_dlpfc))
length(uni_expressed_dlpfc)
# [1] 33899
identical(uni_expressed, uni_expressed_dlpfc)
# [1] FALSE

## Just checking with sacc
load(rse_files['sacc'], verbose = TRUE)
uni_expressed_list_sacc <- list(
    gene = unique(rowRanges(rse_gene)$ensemblID),
    exon = unique(rowRanges(rse_exon)$ensemblID),
    jxn = unique(gsub('\\..*', '', rowRanges(rse_jxn)$newGeneID)),
    tx = unique(gsub('\\..*', '', rowRanges(rse_tx)$gene_id))
)
sapply(uni_expressed_list_sacc, length)
#  gene  exon   jxn    tx
# 25136 23628 19185 27059
identical(sapply(uni_expressed_list, length), sapply(uni_expressed_list_sacc, length))
# [1] TRUE
uni_expressed_sacc <- unique(unlist(uni_expressed_list_sacc))
length(uni_expressed_sacc)
# [1] 33899
identical(uni_expressed, uni_expressed_sacc)
# [1] TRUE

## ok, there are minor differences between amygdala/sacc and DLPFC
uni_expressed <- unique(c(uni_expressed, uni_expressed_dlpfc))
length(uni_expressed)
# [1] 33939

save(uni_expressed, file = 'go_files/uni_expressed.Rdata')


## Load the 881 eQTL tested genes info
eqtl_files <- dir(pattern = 'raggr_suggestive881_snps_.*full.csv')
names(eqtl_files) <- gsub('raggr_suggestive881_snps_|_eqtls_full.csv', '', eqtl_files)
## They are actually the same length, so I only need 1 for the universe
system('wc -l *full.csv')

qtl <- fread(eqtl_files[1])
uni_tested <- unique(qtl$ensemblID)
length(uni_tested)
# [1] 6210

save(uni_tested, file = 'go_files/uni_tested.Rdata')
rm(qtl)

## Load the 881 eQTL results and get the unique gene IDs
eqtl_files <- dir(pattern = 'raggr_881_snps_.*')
names(eqtl_files) <- gsub('raggr_881_snps_|_eqtls_fdr01.csv', '', eqtl_files)

sig_gene <- lapply(eqtl_files, function(ef) {
    message(paste(Sys.time(), 'loading', ef))
    qtl <- fread(ef)
    unique(qtl$ensemblID)
})
sapply(sig_gene, length)
# amyg dlpfc  sacc
#  263   250   381

save(sig_gene, file = 'go_files/sig_gene_881.Rdata')

## Load the 881 eQTL results to subset to FDR < 0.05
eqtl_raw_files <- dir(pattern = 'mergedEqtl.*rda')
names(eqtl_raw_files) <- gsub('.*output_|_raggr.*', '', eqtl_raw_files)

sig_gene_5perc <- lapply(eqtl_raw_files, function(raw_f) {
    message(paste(Sys.time(), 'loading', raw_f))
    load(raw_f, verbose = TRUE)
    
    allEqtl_sub <- subset(allEqtl, FDR < 0.05)
    
    gene_ids <- sapply(unique(allEqtl_sub$Type), function(type) {
        sub_type <- subset(allEqtl_sub, Type == type)
        res <- unique(sub_type$EnsemblGeneID)
        print(paste('Number of eQTL associations at FDR<5%:', nrow(sub_type), '; number of unique gene ids:', length(res)))
        return(res)
    })
    final_res <- unique(unlist(gene_ids))
    print(paste('Number of unique gene ids across all features:', length(final_res), 'for', raw_f))
    return(final_res)
})

# 2018-12-12 11:05:46 loading mergedEqtl_output_amyg_raggr_4features.rda
# Loading objects:
#   allEqtl
# [1] "Number of eQTL associations at FDR<5%: 3724 ; number of unique gene ids: 179"
# [1] "Number of eQTL associations at FDR<5%: 18026 ; number of unique gene ids: 277"
# [1] "Number of eQTL associations at FDR<5%: 10429 ; number of unique gene ids: 233"
# [1] "Number of eQTL associations at FDR<5%: 5051 ; number of unique gene ids: 193"
# [1] "Number of unique gene ids across all features: 645 for mergedEqtl_output_amyg_raggr_4features.rda"
# 2018-12-12 11:06:20 loading mergedEqtl_output_dlpfc_raggr_4features.rda
# Loading objects:
#   allEqtl
# [1] "Number of eQTL associations at FDR<5%: 3470 ; number of unique gene ids: 174"
# [1] "Number of eQTL associations at FDR<5%: 12790 ; number of unique gene ids: 251"
# [1] "Number of eQTL associations at FDR<5%: 6165 ; number of unique gene ids: 189"
# [1] "Number of eQTL associations at FDR<5%: 5589 ; number of unique gene ids: 231"
# [1] "Number of unique gene ids across all features: 613 for mergedEqtl_output_dlpfc_raggr_4features.rda"
# 2018-12-12 11:06:54 loading mergedEqtl_output_sacc_raggr_4features.rda
# Loading objects:
#   allEqtl
# [1] "Number of eQTL associations at FDR<5%: 7364 ; number of unique gene ids: 302"
# [1] "Number of eQTL associations at FDR<5%: 35700 ; number of unique gene ids: 436"
# [1] "Number of eQTL associations at FDR<5%: 13210 ; number of unique gene ids: 303"
# [1] "Number of eQTL associations at FDR<5%: 8610 ; number of unique gene ids: 272"
# [1] "Number of unique gene ids across all features: 912 for mergedEqtl_output_sacc_raggr_4features.rda"

sapply(sig_gene_5perc, length)
# amyg dlpfc  sacc
#  645   613   912
length(unique(unlist(sig_gene_5perc)))
# [1] 1380
save(sig_gene_5perc, file = 'go_files/sig_gene_5perc.Rdata')

run_go <- function(genes_ens, ont = c('BP', 'MF', 'CC'), universe) {
    # ## Change to ENSEMBL ids
    # genes_ens <- sapply(genes, function(x) { gsub('\\..*', '', x) })

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(genes_ens, fun = "enrichGO",
            OrgDb = 'org.Hs.eg.db', universe = universe,
            ont = bp, pAdjustMethod = "BH",
            #pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(universe, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        organism = 'hsa', pAdjustMethod = 'BH', universe = uni_ncbi, 
        #pvalueCutoff = 0.1, qvalueCutoff = 0.05,
        keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}

## FDR < 1%
go_tested <- run_go(sig_gene, universe = uni_tested)
sapply(go_tested, class)
#                     BP                     MF                     CC
# "compareClusterResult"                 "NULL" "compareClusterResult"
#                   KEGG
# "compareClusterResult"
save(go_tested, file = 'go_files/go_tested.Rdata')

go_expressed <- run_go(sig_gene, universe = uni_expressed)
sapply(go_expressed, class)
#                     BP                     MF                     CC
# "compareClusterResult"                 "NULL"                 "NULL"
#                   KEGG
# "compareClusterResult"
save(go_expressed, file = 'go_files/go_expressed.Rdata')


## FDR < 5%
if(FALSE) {
    load('go_files/uni_tested.Rdata', verbose = TRUE)
    load('go_files/uni_expressed.Rdata', verbose = TRUE)
}
go_tested_5perc <- run_go(sig_gene_5perc, universe = uni_tested)
sapply(go_tested_5perc, class)
    # BP     MF     CC
# "NULL" "NULL" "NULL"
save(go_tested_5perc, file = 'go_files/go_tested_5perc.Rdata')

go_expressed_5perc <- run_go(sig_gene_5perc, universe = uni_expressed)
sapply(go_expressed_5perc, class)
#                     BP                     MF                     CC
# "compareClusterResult"                 "NULL" "compareClusterResult"
save(go_expressed_5perc, file = 'go_files/go_expressed_5perc.Rdata')



go_table <- function(go_object) {
    ## Drop any failed results (aka, no enrichment)
    go_object <- go_object[sapply(go_object, class) == 'compareClusterResult']
    if(length(go_object) == 0) return(NULL)

    message(paste(Sys.time(), 'found enrichment for the ontologies', paste(names(go_object), collapse = ', ') ))

    ## Convert to a table (aka, flatten the results)
    go_table <- do.call(rbind, mapply(function(x, name_x) {
        result <- as.data.frame(x)
        result$ontology <- name_x
        return(result)
    }, go_object, names(go_object), SIMPLIFY = FALSE))

    ## Then to a DataFrame
    go_table_df <- DataFrame(go_table)
    go_table_df$geneID <- CharacterList(strsplit(go_table$geneID, '/'))
    
    return(go_table_df)
}

## From https://github.com/LieberInstitute/brainseq_phase2/blob/ddc659bcdda712ae653b8f744299438b7483315c/supp_tabs/create_supp_tables.R#L110-L118
fix_csv <- function(df) {
    for(i in seq_len(ncol(df))) {
        if(any(grepl(',', df[, i]))) {
            message(paste(Sys.time(), 'fixing column', colnames(df)[i]))
            df[, i] <- gsub(',', ';', df[, i])
        }
    }
    return(df)
}

go_write <- function(go_table, go_file) {
    go_df <- fix_csv(as.data.frame(go_table))
    write.csv(go_df, go_file, quote = FALSE)
}

## FDR < 1%
go_tested_table <- go_table(go_tested)
# 2018-12-11 15:41:35 found enrichment for the ontologies BP, CC, KEGG
save(go_tested_table, file = 'go_files/go_tested_table.Rdata')
go_write(go_tested_table, 'go_files/go_tested_table.csv')
# 2018-12-11 15:42:04 fixing column geneID

go_expressed_table <- go_table(go_expressed)
# 2018-12-11 16:01:28 found enrichment for the ontologies BP, KEGG
save(go_expressed_table, file = 'go_files/go_expressed_table.Rdata')
go_write(go_expressed_table, 'go_files/go_expressed_table.csv')
# 2018-12-11 16:01:30 fixing column geneID

## FDR < 5%
go_tested_5perc_table <- go_table(go_tested_5perc)
save(go_tested_5perc_table, file = 'go_files/go_tested_5perc_table.Rdata')
go_write(go_tested_5perc_table, 'go_files/go_tested_5perc_table.csv')

go_expressed_5perc_table <- go_table(go_expressed_5perc)
# 2018-12-12 12:08:00 found enrichment for the ontologies BP, CC
save(go_expressed_5perc_table, file = 'go_files/go_expressed_5perc_table.Rdata')
go_write(go_expressed_5perc_table, 'go_files/go_expressed_5perc_table.csv')
# 2018-12-12 12:08:10 fixing column geneID


## Make the plots: requires re-loading using:
# module load conda_R/3.4.x
load('go_files/go_tested.Rdata', verbose = TRUE)
load('go_files/go_expressed.Rdata', verbose = TRUE)
load('go_files/go_tested_5perc.Rdata', verbose = TRUE)
load('go_files/go_expressed_5perc.Rdata', verbose = TRUE)

simplify_go <- function(x) {
    # gsub('amygdala', 'amyg', gsub('dep', 'de', gsub('ptsd', 'pt', x)))
    ## Don't really need this function for this project
    x
}

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        go@compareClusterResult$Cluster <- simplify_go(go@compareClusterResult$Cluster)
        names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}

## FDR < 1%
pdf('go_files/go_tested.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_tested)
dev.off()

pdf('go_files/go_tested_all.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_tested, cat = NULL)
dev.off()

pdf('go_files/go_expressed.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_expressed)
dev.off()

pdf('go_files/go_expressed_all.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_expressed, cat = NULL)
dev.off()

## FDR < 5%
pdf('go_files/go_tested_5perc.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_tested_5perc)
dev.off()

pdf('go_files/go_tested_5perc_all.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_tested_5perc, cat = NULL)
dev.off()

pdf('go_files/go_expressed_5perc.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_expressed_5perc)
dev.off()

pdf('go_files/go_expressed_5perc_all.pdf', width = 9, height = 9, useDingbats = FALSE)
plot_go(go_expressed_5perc, cat = NULL)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2018-12-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  AnnotationDbi        * 1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.2    2018-11-28 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  clusterProfiler      * 3.10.0    2018-10-30 [1] Bioconductor
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 [2] CRAN (R 3.5.0)
#  cowplot                0.9.3     2018-07-15 [1] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table           * 1.11.8    2018-09-30 [1] CRAN (R 3.5.1)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  DO.db                  2.9       2018-05-03 [1] Bioconductor
#  DOSE                   3.8.0     2018-10-30 [1] Bioconductor
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  enrichplot             1.2.0     2018-10-30 [2] Bioconductor
#  europepmc              0.3       2018-04-20 [2] CRAN (R 3.5.1)
#  farver                 1.1.0     2018-11-20 [1] CRAN (R 3.5.1)
#  fastmatch              1.1-0     2017-01-28 [1] CRAN (R 3.5.0)
#  fgsea                  1.8.0     2018-10-30 [1] Bioconductor
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggforce                0.1.3     2018-07-07 [2] CRAN (R 3.5.1)
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  ggplotify              0.0.3     2018-08-03 [2] CRAN (R 3.5.1)
#  ggraph                 1.0.2     2018-07-07 [2] CRAN (R 3.5.1)
#  ggrepel                0.8.0     2018-05-09 [1] CRAN (R 3.5.0)
#  ggridges               0.5.1     2018-09-27 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  GO.db                  3.7.0     2018-11-02 [1] Bioconductor
#  GOSemSim               2.8.0     2018-10-30 [1] Bioconductor
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 3.5.0)
#  gridGraphics           0.3-0     2018-04-03 [2] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5     2018-07-19 [2] CRAN (R 3.5.1)
#  httr                   1.3.1     2017-08-20 [1] CRAN (R 3.5.0)
#  igraph                 1.2.2     2018-07-27 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jsonlite               1.5       2017-06-01 [2] CRAN (R 3.5.0)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  MASS                   7.3-51.1  2018-11-01 [3] CRAN (R 3.5.1)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  org.Hs.eg.db         * 3.7.0     2018-11-02 [1] Bioconductor
#  pillar                 1.3.0     2018-07-14 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  qvalue                 2.14.0    2018-10-30 [1] Bioconductor
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.0.1   2018-10-25 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  rvcheck                0.1.3     2018-12-06 [1] CRAN (R 3.5.1)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 1.4.2     2018-01-22 [1] CRAN (R 3.5.0)
#  tidyr                  0.8.2     2018-10-28 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  triebeard              0.3.0     2016-08-04 [1] CRAN (R 3.5.0)
#  tweenr                 1.0.0     2018-09-27 [1] CRAN (R 3.5.1)
#  units                  0.6-2     2018-12-05 [1] CRAN (R 3.5.1)
#  UpSetR                 1.3.3     2017-03-21 [1] CRAN (R 3.5.0)
#  urltools               1.7.1     2018-08-03 [1] CRAN (R 3.5.1)
#  viridis                0.5.1     2018-03-29 [2] CRAN (R 3.5.0)
#  viridisLite            0.3.0     2018-02-01 [2] CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  xml2                   1.2.0     2018-01-24 [2] CRAN (R 3.5.0)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
# >
