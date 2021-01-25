library('recount.bwtool')
library('SummarizedExperiment')
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'region', 'r', 1, 'character', 'Which brain region to subset',
    'casestatus', 'c', 1, 'character', 'Dx status. Either case, control or all',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list('region' = 'amygdala', 'casestatus' = 'all')
    opt <- list('region' = 'sacc', 'casestatus' = 'all')
    opt <- list('region' = 'dlpfc', 'casestatus' = 'all')
}

stopifnot(opt$region %in% c('amygdala', 'sacc', 'dlpfc'))
stopifnot(opt$casestatus %in% c('case', 'control', 'all'))

## Are we on JHPCE? Print some helpful info
jhpce <- grepl('compute-', Sys.info()['nodename'])
if(jhpce) {
    message(paste(Sys.time(), 'Remember to run before this R script:
    module load ucsctools
    module load wiggletools/default
    '))
}

## Load pheno data
f <- paste0('/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_', opt$region, '.rda')
stopifnot(file.exists(f))
message(paste(Sys.time(), 'loading', f))
load(f, verbose = TRUE)

## Remove unwanted stuff
rm(rse_exon, rse_jxn, rse_tx)

## Set variables that change by region
if(opt$region %in% c('amygdala', 'sacc')) {
    bwroot <- '/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/Coverage/'
    strands <- c('.Forward', '.Reverse')
    ids <- colData(rse_gene)$SAMPLE_ID
    dx <- as.character(colData(rse_gene)$PrimaryDx)
    dx <- ifelse(dx %in% c('Bipolar', 'Other'), 'Bipolar', ifelse(dx == 'Control', 'Control', NA))
} else if (opt$region == 'dlpfc') {
    bwroot <- '/dcl01/ajaffe/data/lab/brainseq_phase1/preprocessed_data/Coverage/'
    dx <- 'Dx'
    strands <- c('')
    ids <- sapply(colData(rse_gene)$SAMPLE_ID, '[', 1)
    dx <- colData(rse_gene)$Dx
}

## Subset to samples of interest
if(opt$casestatus != 'all') {
    if(opt$casestatus == 'case') {
        ids <- ids[dx == 'Bipolar']
    } else {
        ids <- ids[dx == 'Control']
    }
}

## Find BigWig files and compute mean by strand
sapply(strands, function(strand) {
    message(paste(Sys.time(), 'processing strand', strand))
    bws <- paste0(bwroot, ids, strand, '.bw')
    bws_exist <- file.exists(bws)
    message(paste(Sys.time(), 'found', sum(bws_exist), 'BigWig files'))
    stopifnot(all(bws_exist))
    
    outfile <- paste0('mean_', opt$region, '_', opt$casestatus, '_n', length(bws), strand)
    tmp_dir <- paste0('temp_', opt$region, '_', opt$casestatus, strand)
    dir.create(tmp_dir, showWarnings = FALSE)
    compute_mean(bws, outfile = outfile, tempdir = tmp_dir)
})

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
