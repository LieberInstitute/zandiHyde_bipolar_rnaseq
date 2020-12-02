library(tidyverse)
library(data.table)
library(SummarizedExperiment)

data.table::setDTthreads(threads = 1)

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
    rse
)

rowData(sacc_rse)$EntrezID
rowData(amygdala_rse)$EntrezID

