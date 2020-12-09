library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(devtools)

data.table::setDTthreads(threads = 1)

## Define background universe of genes ####

load("rda/generate_plots_data.RData")

load("amygdala_gene/working_rse.RData",
     verbose = TRUE)

amygdala_rse <- rse

load("sacc_gene/working_rse.RData",
     verbose = TRUE)

sacc_rse <- rse

rse <- NULL

# remove extra objects
rm(
    amyg_rho,
    axisdf,
    don,
    don_key,
    fin_plot,
    i,
    intctv_plot,
    p,
    sacc_rho,
    statOutGene,
    twas_exp_fin,
    twas_var,
    twas_z_wide,
    twas_z_sig_tables,
    rse,
    twas_z
)

## Find Entrez IDs ####
# create a new EntrezID column in twas_z_sacc which contains the Entrez IDs of sACC in the rse
twas_z_sacc[, EntrezID := rowData(sacc_rse)[match(twas_z_sacc$geneid, rowData(sacc_rse)$gencodeID),]$EntrezID]

twas_z_amyg[, EntrezID := rowData(amygdala_rse)[match(twas_z_amyg$geneid, rowData(amygdala_rse)$gencodeID),]$EntrezID]

# > summary(is.na(twas_z_sacc$EntrezID))
#    Mode   FALSE    TRUE
# logical    7246    2754

# > summary(is.na(twas_z_amyg$EntrezID))
#    Mode   FALSE    TRUE
# logical    6437    2538

## Remove NAs ####

twas_z_amyg <- twas_z_amyg[!is.na(twas_z_amyg$EntrezID),]
twas_z_sacc <- twas_z_sacc[!is.na(twas_z_sacc$EntrezID),]

twas_z_both <- rbind(twas_z_amyg, twas_z_sacc)
twas_z_both <- twas_z_both[!is.na(twas_z_both$EntrezID),]

## Calculate FDR < 5% ####

twas_z_both$fdr.p <- p.adjust(twas_z_both$TWAS.P, 'fdr')
twas_z_amyg$fdr.p <- p.adjust(twas_z_amyg$TWAS.P, 'fdr')
twas_z_sacc$fdr.p <- p.adjust(twas_z_sacc$TWAS.P, 'fdr')

twas_z_both_fdr <- twas_z_both[fdr.p < 0.05,]
twas_z_amyg_fdr <- twas_z_amyg[fdr.p < 0.05,]
twas_z_sacc_fdr <- twas_z_sacc[fdr.p < 0.05,]

## clusterProfiler ####

# compareCluster

# test_twas_z_both_fdr <- twas_z_both_fdr[1:100]
# test_twas_z_both <- twas_z_both[1:100]
#
# go_test_both <- enrichGO(gene = test_twas_z_both_fdr$EntrezID, OrgDb = "org.Hs.eg.db",
#          keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 1,
#          pAdjustMethod = "fdr", universe = test_twas_z_both$EntrezID,
#          qvalueCutoff = 1, readable = TRUE)

go_both <- enrichGO(gene = twas_z_both_fdr$EntrezID, OrgDb = "org.Hs.eg.db",
         keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 1,
         pAdjustMethod = "fdr", universe = twas_z_both$EntrezID,
         qvalueCutoff = 1, readable = TRUE)

go_amyg <- enrichGO(gene = twas_z_both_fdr$EntrezID, OrgDb = "org.Hs.eg.db",
         keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 1,
         pAdjustMethod = "fdr", universe = twas_z_both$EntrezID,
         qvalueCutoff = 1, readable = TRUE)

go_sacc <- enrichGO(gene = twas_z_both_fdr$EntrezID, OrgDb = "org.Hs.eg.db",
         keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 1,
         pAdjustMethod = "fdr", universe = twas_z_both$EntrezID,
         qvalueCutoff = 1, readable = TRUE)

save(go_both, file = "rda/go_both_enrichment.rda")
save(go_amyg, file = "rda/go_amyg_enrichment.rda")
save(go_sacc, file = "rda/go_sacc_enrichment.rda")

session_info()
