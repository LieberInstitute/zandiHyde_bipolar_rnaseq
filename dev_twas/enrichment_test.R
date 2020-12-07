library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(clusterProfiler)

data.table::setDTthreads(threads = 1)

## Define background universe of genes ####

load("generate_plots_data.RData")

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

go_both <- compareCluster(twas_z_both_fdr$fdr.p, univ = twas_z_both$fdr.p,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)

go_amyg <- compareCluster(twas_z_amyg_fdr$fdr.p, univ = twas_z_amyg$fdr.p,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)

go_sacc <- compareCluster(twas_z_sacc_fdr$fdr.p, univ = twas_z_sacc$fdr.p,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)

save(go_both, file = "go_both_enrichment.rda")
save(go_amyg, file = "go_amyg_enrichment.rda")
save(go_sacc, file = "go_sacc_enrichment.rda")
