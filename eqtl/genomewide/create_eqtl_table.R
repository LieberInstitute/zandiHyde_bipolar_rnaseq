##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("mergedEqtl_output_sacc_genomewide_4features_FDR01.rda", verbose=TRUE)
sacc = allEqtlFDR01
load("mergedEqtl_output_amyg_genomewide_4features_FDR01.rda", verbose=TRUE)
amyg = allEqtlFDR01
load("mergedEqtl_output_dlpfc_genomewide_4features_FDR01.rda", verbose=TRUE)
dlp = allEqtlFDR01

################
## metrics

## total features
nrow(sacc)  ## 12,345,494
nrow(amyg)   ## 9,292,449
nrow(dlp) ## 8,331,329

## per feature
table(sacc$Type)
#    Exon    Gene      Jxn       Tx
# 6,560,649  1,241,720 3,087,320 1,455,805
table(amyg$Type)
#      Exon     Gene        Jxn      Tx
# 4,866,350  937,203 2,341,997 1,146,899
table(dlp$Type)
#     Exon      Gene       Jxn      Tx
# 4,322,899  780,312 2,112,521 1,115,597

## unique ensemblIDs
tapply(sacc$EnsemblGeneID, sacc$Type, function(x) length(unique(x)))
# Exon    Gene  Jxn   Tx
#  12115  9904  9564  8489
tapply(amyg$EnsemblGeneID, amyg$Type, function(x) length(unique(x)))
# Exon    Gene  Jxn   Tx
#  10507  8034  8421  7347
tapply(dlp$EnsemblGeneID, dlp$Type, function(x) length(unique(x)))
# Exon  Gene  Jxn   Tx
#  9852 7487 8198 7639


	
	
## unique SNPs
length(unique(sacc$snps))  ## 1,321,958
length(unique(amyg$snps))  ## 1,096,532
length(unique(dlp$snps))   ##  985,686



length(unique(sacc$gene))  ##  126,046
length(unique(amyg$gene))  ##  99,218
length(unique(dlp$gene))   ##  94,934

## unique features
tapply(sacc$gene, sacc$Type, function(x) length(unique(x)))
#  Exon  Gene   Jxn    Tx
# 70000  9904 32478 13664
tapply(amyg$gene, amyg$Type, function(x) length(unique(x)))
#  Exon  Gene   Jxn    Tx
# 54301  8034 25499 11384
tapply(dlp$gene, dlp$Type, function(x) length(unique(x)))
#  Exon  Gene   Jxn    Tx
# 50927  7487 24639 11881


## SNP-feature pairs
nrow(sacc)  ## 12,345,494
nrow(amyg)   ## 9,292,449
nrow(dlp) ## 8,331,329
table(sacc$Type)
table(amyg$Type)
table(dlp$Type)


## Unique symbols in SNP-feature pairs
length(unique(sacc$Symbol))  ## 15705
length(unique(amyg$Symbol))  ## 14196
length(unique(dlp$Symbol))  ## 15092
tapply(sacc$Symbol, sacc$Type, function(x) length(unique(x)))
tapply(amyg$Symbol, amyg$Type, function(x) length(unique(x)))
tapply(dlp$Symbol, dlp$Type, function(x) length(unique(x)))

# > tapply(sacc$Symbol, sacc$Type, function(x) length(unique(x)))
 # Exon  Gene   Jxn    Tx
# 10526  7927  9387  8443
# > tapply(amyg$Symbol, amyg$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 9142 6331 8263 7305
# > tapply(dlp$Symbol, dlp$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 9825 7454 8190 7590

	
	
	
	
	


