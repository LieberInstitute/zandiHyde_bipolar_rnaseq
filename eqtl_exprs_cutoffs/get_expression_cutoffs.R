####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(segmented)

## load
load("../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseGene_n511.rda")

## recount getRPKM version ##
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
getRPM = function(rse, target = 80e6) {
	require(SummarizedExperiment)
	mapped <- colSums(assays(rse)$counts) 
	bg = matrix(rep(mapped/target), nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

####################
### Regions combined
exprs <- list(
    'Gene' = getRPKM(rse_gene, length_var = "Length"),
    'Exon' = getRPKM(rse_exon, length_var = "Length"),
    'Jxn' = getRPM(rse_jxn, target=10e6),
    'Tx' = assays(rse_tx)$tpm
)
## Identify potential cutoffs
seed <- 20171026
seeds <- seed + 0:3
names(seeds) <- names(exprs)
cutoffs <- sapply(names(exprs), function(type) {  
    message(type)
    # pdf(paste0('suggested_expr_cutoffs_', tolower(type), '.pdf'), width = 12)
    cuts <- expression_cutoff(exprs[[type]], seed = seeds[type])
    message(paste(cuts, collapse = ' '))
    cut <- max(cuts)
    # dev.off()   
    return(cut)
})
# Gene
# 2017-12-21 14:26:53 the suggested expression cutoff is 0.26 0.17
# Exon
# 2017-12-21 14:31:22 the suggested expression cutoff is 0.3 0.21
# Jxn
# 2017-12-21 14:35:11 the suggested expression cutoff is 0.21 0.33
# Tx
# 2017-12-21 14:37:48 the suggested expression cutoff is 0.4 0.24


### Final cutoffs used:
# Gene 0.25
# Exon 0.30
# Jxn 0.35
# Tx 0.40










