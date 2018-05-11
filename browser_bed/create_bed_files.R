## Based on:
# /users/ajaffe/Lieber/Projects/RNAseq/Astrazeneca/Round3/round3_bedTrack.R

library('rtracklayer')
library('RColorBrewer')
library('jaffelab')
library('devtools')

## Load eQTL results from the 31 main SNPs and their proxy SNPs
efiles <- dir('../eqtl/raggr_gwSignificant', pattern = 'raggr_31_snps_', full.names = TRUE)
names(efiles) <- ss(gsub('_eqtls_fdr01.csv', '', efiles), '_snps_', 2)
etabs <- lapply(efiles, read.csv, header = TRUE)

get_col <- function(pval) {
    bedShading <- cut(pval,
    	breaks = c(0,3, 5, 8, 10, 12, 15, 20, 1000), label=FALSE,
    	include.lowest=TRUE)
        
    Ngroup = max(bedShading)
    pal = brewer.pal(7,"RdBu")
    posCols = colorRampPalette(pal[4:7])(Ngroup)
    names(posCols) = 1:Ngroup
    
    negCols = colorRampPalette(pal[1:4])(Ngroup)
    names(negCols) = seq(-1*Ngroup,-1)
    cols = c(posCols, negCols)
    
    bedColors = cols[match(as.character(bedShading), names(cols))]
    
    tmp = t(col2rgb(bedColors))
    paste(tmp[,1], tmp[,2], tmp[,3], sep= ",")
}

ebed <- lapply(etabs, function(x) {
    df <- data.frame(
        'chr' = x$chr_hg38,
        'start' = x$feat_start,
        'end' = x$feat_end,
        'name' = x$SNP,
        'score' = -log10(x$FDR),
        'strand' = x$strand,
        'thickStart' = x$feat_start,
        'thickEnd' = x$feat_end,
        'itemRgb' = get_col(-log10(x$FDR)),
        'status' = x$Status,
        'type' = x$Type,
        stringsAsFactors = FALSE
    )
    
    df[order(df$score, decreasing = TRUE), ]
    
    split(df[, -which(colnames(df) %in% c('status', 'type'))], paste0(df$status, '_', df$type))
})

stopifnot(all(sapply(ebed, length) == 8))

get_header <- function(region, typestatus) {
    paste0("track name=ZandiBipolar_eQTL_", region, '_', typestatus,
    		" description='ZandiBipolar eQTL hits - ", region, ', ',
            ss(typestatus, '_'), ', ', ss(typestatus, '_', 2),
    		"' visibility=2 itemRgb='On'")
}

## Write to files
dir.create('bed', showWarnings = FALSE)
xx <- lapply(names(ebed), function(reg) {
    lapply(names(ebed[[1]]), function(typstat) {
        bed <- paste0('bed/zandiBipolar_31main_bedTrack_05112018_', reg, '_', typstat, '.bed')
        
    	write.table(get_header(reg, typstat), file = bed,
    		row.names=FALSE, col.names=FALSE, quote=FALSE)
        xx <- ebed[[reg]][[typstat]]
        if(any(xx$strand == '*')) xx$strand[xx$strand == '*'] <- '+'
    	write.table(xx, file = bed,
    		row.names=FALSE, col.names=FALSE, quote=FALSE,append=TRUE)
        return(NULL)
    })
})

## Compress into a single tar ball for sharing
system('tar -zcvf bed.tar.gz bed')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
