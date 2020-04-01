##
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")

### make "Other" dx bipolar
pd = colData(rse_gene)
pd$PrimaryDx[pd$PrimaryDx=="Other"] = "Bipolar"
pd$PrimaryDx = factor(pd$PrimaryDx, levels = c("Control", "Bipolar"))

## load SNP data
load("../../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum



# statistical model ##
load("../rdas/pcs_4features_amyg.rda")
	
modG = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + genePCs, data = pd)
modE = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + exonPCs, data = pd)
modJ = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + jxnPCs, data = pd)
modT = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]) + txPCs, data = pd)

################
## residualize expression	
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
cleanGeneSub = log2(geneRpkm+1)		
cleanGeneSub = cleaningY(cleanGeneSub, modG, P=1)

cleanExonSub = log2(exonRpkm+1)		
cleanExonSub = cleaningY(cleanExonSub, modE, P=1)

cleanJxnSub = log2(jxnRp10m+1)		
cleanJxnSub = cleaningY(cleanJxnSub, modJ, P=1)

cleanTxSub = log2(txTpm+1)		
cleanTxSub = cleaningY(cleanTxSub, modT, P=1)


######################
##### Conditional ####
######################

## pgc rank
pgc = read.csv("../../../PGC_risk_loci.csv")
pgc = pgc[order(pgc$P.value),]
pgc$Rank = 1:nrow(pgc)

## eqtl results
pgcEqtlsFinal = read.csv("raggr_881_snps_amyg_eqtls_fdr01.csv")
names(pgcEqtlsFinal)[c(2,6)] = c("SNP","Feature")
pgcEqtlsFinal$pgcFinalRank = pgc$Rank[match(pgcEqtlsFinal$IndexSNP_hg19POS,pgc$hg19POS)]

pgcLead = pgcEqtlsFinal[which(pgcEqtlsFinal$leadVariant_indicator == TRUE),]
pgcLead = pgcLead[order(pgcLead$pvalue),]

## filter expression
yExprs = rbind(cleanGeneSub[rownames(cleanGeneSub) %in% pgcEqtlsFinal$Feature,],
	cleanExonSub[rownames(cleanExonSub) %in% pgcEqtlsFinal$Feature,],
	cleanJxnSub[rownames(cleanJxnSub) %in% pgcEqtlsFinal$Feature,],
	cleanTxSub[rownames(cleanTxSub) %in% pgcEqtlsFinal$Feature,])
yExprs = yExprs[match(pgcEqtlsFinal$Feature, rownames(yExprs)),]

s = snp[match(pgcEqtlsFinal$SNP, rownames(snp)),]

##### analyses by region, post hoc on clean to compare
statPost = t(sapply(1:nrow(pgcEqtlsFinal), function(i) {
	summary(lm(yExprs[i,] ~ s[i,]))$coef[2,-2]
}))
statPost = as.data.frame(statPost)
colnames(statPost) = c("beta", "statistic", "pvalue")

plot(statPost$statistic ~ pgcEqtlsFinal$statistic) # looks okay

rm(list=ls()[! ls() %in% c("pgcEqtlsFinal","pgcLead","s","yExprs")])

###### split by region
rIndexes = splitit(pgcLead$pgcFinalRank)
sigIndexes = sapply(rIndexes, function(ii) {
	cat(".")
	best = ii[1]
	rest = ii[-1]
	theSnp = s[best,]
	while(length(rest) > 1) { 
		if(length(best) == 1) z = yExprs[best,] else z = t(yExprs[best,])
		coefs = as.data.frame(t(sapply(rest, function(i) summary(lm(yExprs[i,] ~ 
			theSnp + z))$coef[2,-2])))
		colnames(coefs) = c("beta", "statistic", "pvalue")
		
		if(sum(coefs$pvalue < 0.05) == 0) {
			return(best) 
		} else {
			best = c(best, rest[which(coefs$pvalue < 0.05)[1]])
			rest = rest[which(coefs$pvalue < 0.05)[-1]]
		}
	}
	
	return(best)
})
indeptIndex = unlist(sigIndexes)

pgcLead$condIndept = 0
pgcLead$condIndept[indeptIndex] = 1


########################
## check results out ###
########################

eqtl = pgcLead

### bring in case-control DE
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint.rda")

### put in order
statOut = statOut[match(eqtl$Feature,rownames(statOut)),]

### add dir and pval from case-control DE
statOut$bipolar_feat_dir = ifelse(statOut$logFC_Amyg > 0, "UP", "DOWN")
eqtl$bipolar_feat_dir = statOut$bipolar_feat_dir
eqtl$bipolar_feat_pval = statOut$P.Value_Amyg

eqtl = eqtl[order(eqtl$pgcFinalRank, eqtl$pvalue),]

# write.csv(eqtl, file="pgc_conditional_by_leadVar_amyg.csv")





####################################
#### transcript plots
library(GenomicFeatures)
library('BSgenome.Hsapiens.UCSC.hg38')
library(derfinderPlot)

## GENCODE
chrInfo = getChromInfoFromUCSC("hg38")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

gencode_v25 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
seqlevels(gencode_v25, pruning.mode="coarse") = paste0("chr", c(1:22,"X","Y"))
seqinfo(gencode_v25) = si
gencode_v25_txdb = makeTxDbFromGRanges(gencode_v25)
txByExon = exonsBy(gencode_v25_txdb, use.names=TRUE)

## select gene coords
# chrom = geneInfo$seqnames
# minRange = geneInfo$start
# maxRange = geneInfo$end
# strnd=geneInfo$strand
# gr = GRanges(chrom, IRanges(minRange,maxRange), strnd)
gr = GRanges("chr11", IRanges(66159000,66244750))

pdf("test_PACS1_structure.pdf",h=3.5,w=7)
par(mar=c(5,3,2,2))
ov <- findOverlaps(gr, txByExon)
txList <- split(txByExon[subjectHits(ov)], queryHits(ov))
poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
yrange <- range(poly.data$exon$y)
xrange <- c(start(gr),end(gr))
plot(0,0, type="n", xlim = xrange , ylim = yrange + c(-0.75, 0.75),
	yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
	cex.axis = 1.3, cex.lab =1.3, cex.main =1.3, main="PACS1")
rect(xleft=66220288, xright=66220791, ybottom=-100, ytop=100, border=FALSE, col="peachpuff")
box(which = "plot") ## fix border
derfinderPlot:::.plotPolygon(poly.data$exon, 'blue')
derfinderPlot:::.plotPolygon(poly.data$intron, 'lightblue')
abline(v=66177715, lty=2, col="red")  ## SNP
yTx <- unique(yrange / abs(yrange))
if(length(yTx) > 1) {
	axis(2, yTx, c('-', '+')[c(-1, 1) %in% yTx], 
		tick = FALSE, las = 1, cex.axis = 3)
	abline(h = 0, lty = 3)
}
dev.off()
	


### hg19:
# library(rtracklayer)
# library("EnsDb.Hsapiens.v75")
# ensTxDb = EnsDb.Hsapiens.v75
# seqlevelsStyle(ensTxDb) <- "UCSC"
	
# ensTxDb = loadDb("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/ensembl_v75_txdb_update.sqlite")
# ensTxDb = keepStandardChromosomes(ensTxDb, species = "Homo_sapiens")
# ensTxDb = renameSeqlevels(ensTxDb, paste0("chr", seqlevels(ensTxDb)))
# txByExon <- exonsBy(ensTxDb)

# ## select gene coords
# gr = GRanges("chr19", IRanges(19387300,19431350))

# pdf("test_SUGP1_structure.pdf",h=3.5,w=7)
# par(mar=c(5,2,2,2))
# ov <- findOverlaps(gr, txByExon)
# txList <- split(txByExon[subjectHits(ov)], queryHits(ov))
# poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
# yrange <- range(poly.data$exon$y)
# xrange <- c(start(gr),end(gr))
# plot(0,0, type="n", xlim = xrange , ylim = yrange + c(-0.75, 0.75),
	# yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
	# cex.axis = 1.3, cex.lab =1.3, cex.main =1.3, main="SUGP1")
# derfinderPlot:::.plotPolygon(poly.data$exon, 'blue')
# derfinderPlot:::.plotPolygon(poly.data$intron, 'lightblue')
# yTx <- unique(yrange / abs(yrange))
# if(length(yTx) > 1) {
	# axis(2, yTx, c('-', '+')[c(-1, 1) %in% yTx], 
		# tick = FALSE, las = 1, cex.axis = 3)
	# abline(h = 0, lty = 3)
# }
# dev.off()
	
	







########################
## check results out ###
########################


# ### bring in case-control
# load("/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/all_de_features.rda")
# outStatsSub = unlist(endoapply(outStats, function(x) x[names(x) %in% pgcEqtlsFinal$Feature]))
# outStatsSub$Type = ss(names(outStatsSub), "\\.",1)
# outStatsSub$Feature = ss(names(outStatsSub), "\\.",2)
# outStatsMatch = outStatsSub[match(pgcEqtlsFinal$Feature, outStatsSub$Feature)]

# pgcEqtlsFinal$SzDir_qSVA = ifelse(outStatsMatch$log2FC_qsva > 0, "UP", "DOWN")
# pgcEqtlsFinal$SzPval_qSVA = outStatsMatch$pval_qsva

# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA, pgcEqtlsFinal$riskDir))
# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$GTEx_MetaPval < 1e-8], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$GTEx_MetaPval < 1e-8]))
# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$bonf < 0.05], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$bonf < 0.05]))

# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$SzPval_qSVA < 0.1], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$SzPval_qSVA < 0.1]))	
	
# ### clean up
# pgcEqtlsFinal$WhichTx[pgcEqtlsFinal$Type == "Gene"] = NA
# rownames(pgcEqtlsFinal) = NULL
# pgcEqtlsFinal$flipSlopes = NULL
# pgcEqtlsFinal$pgcRow = NULL

# ## write table
# pgcOut = as.data.frame(pgcEqtlsFinal)
# pgcOut$WhichTx = sapply(pgcOut$WhichTx, paste, collapse=";")
# pgcOut = pgcOut[order(pgcOut$pgcFinalRank, pgcOut$pvalue),]
# write.csv(pgcOut, file = "tables/pgcEqtls_fdr01_replStats_withIndept.csv",
	# row.names=FALSE)

# save(pgcEqtlsFinal, file = "rdas/PGC_SZ_indexSnps_annotated_conditionalEqtls.rda")





# #### metrics
# load("rdas/PGC_SZ_indexSnps_annotated_conditionalEqtls.rda")

pgcEqtlsIndept = eqtl[eqtl$condIndept == 1,] 

# ## fix as3mt
# pgcEqtlsIndept$EnsemblID[pgcEqtlsIndept$Symbol == "C10orf32-ASMT"] = "ENSG00000270316"
# pgcEqtlsIndept$EnsemblID[pgcEqtlsIndept$Symbol == "AS3MT-C10orf32-ASMT"] = "ENSG00000270316"
# pgcEqtlsIndept$Symbol[pgcEqtlsIndept$Symbol == "C10orf32-ASMT"] = "AS3MT"
# pgcEqtlsIndept$Symbol[pgcEqtlsIndept$Symbol == "AS3MT-C10orf32-ASMT"] = "AS3MT"

# tableCompare = sapply(list(FDR = eqtl, Indept= pgcEqtlsIndept), function(x) {
	# ttNovel = table(x$pgcFinalRank, x$Class)	## novelty
	# ttType = table(x$pgcFinalRank, x$Type)	## type

	# data.frame(numGwas = length(unique(x$pgcFinalRank)),
		# numGwasNoGene = sum(ttType[,"Gene"] == 0),
		# numGwasNoTxOrGene = sum(rowSums(ttType[,c("Gene", "Tx")]) == 0),
		# numGwasUn = sum(rowSums(ttNovel[,-3]) > 0),
		# numGwasOnlyUn = sum(ttNovel[,"InGen"] == 0),
		# numGwasTx = sum(lengths(sapply(split(x$WhichTx, x$pgcFinalRank), 
			# function(y) unlist(unique(y)))) <= 1))
# })
# tableCompare

## summary stats
dim(pgcEqtlsIndept)
table(pgcEqtlsIndept$Type)
length(unique(pgcEqtlsIndept$ensemblID[!grepl("-", pgcEqtlsIndept$ensemblID) & 
	!is.na(pgcEqtlsIndept$ensemblID) ]))
length(unique(pgcEqtlsIndept$Symbol[!grepl("-", pgcEqtlsIndept$Symbol) & 
	!is.na(pgcEqtlsIndept$Symbol) ]))
	
# number genes
table(lengths(sapply(split(eqtl, eqtl$pgcFinalRank), function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID) & !grepl("-", x$ensemblID)]))))
table(lengths(sapply(split(pgcEqtlsIndept, pgcEqtlsIndept$pgcFinalRank), function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID) & !grepl("-", x$ensemblID)]))))


locusList = split(pgcEqtlsIndept, pgcEqtlsIndept$pgcFinalRank)

## number of genes
table(lengths(sapply(locusList, function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID)]))))


# num tx
t(sapply(locusList, function(x) 
	c(length(unique(x$ensemblID[!is.na(x$ensemblID)])), 
		length(unique(unlist(x$WhichTx[!is.na(x$WhichTx) & x$WhichTx!=""]))))))
		
checkInd0 = which(lengths(sapply(locusList, function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID)]))) == 0)
as.list(locusList[names(checkInd0)])

checkInd1 = which(lengths(sapply(locusList, function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID)]))) == 1)
as.list(locusList[names(checkInd1)])

checkInd2 = which(lengths(sapply(locusList, function(x) 
	unique(x$ensemblID[!is.na(x$ensemblID)]))) == 2)
as.list(locusList[names(checkInd2)])

## one and two gene loci
oneOrTwo = sapply(locusList, function(x) 
	length(unique(x$EnsemblID)) %in% 1:2)
tmp = sapply(locusList[oneOrTwo], function(x) unique(x$EnsemblID))
outTable = data.frame(locus = rep(names(tmp), lengths(tmp)), EnsemblID = unlist(tmp),
	stringsAsFactors=FALSE)
outTable$Gene = pgcEqtlsIndept$Symbol[match(outTable$EnsemblID, pgcEqtlsIndept$EnsemblID)]
outTable$Gene[is.na(outTable$EnsemblID)] = "Intergenic"
outTable$EnsemblID[is.na(outTable$EnsemblID)] = "Intergenic"
outTable$Gene[outTable$Gene == ""] = c("ZSCAN26", "AC068831.1", "AC011816.1", "AC117382.2",
	"SNORA77", "AC090568.2", "AP003049.1")
outTable$SNP = pgcEqtlsIndept$snpRsNum[match(outTable$locus, pgcEqtlsIndept$pgcFinalRank)]
outTable = outTable[,c(1,4,3,2)]

write.csv(outTable, file="tables/table2_eqtl_pgc_oneOrTwo.csv",row.names=FALSE)
	
###############
# make plots ##
###############



	
	
	
# TranscriptDb = makeTxDbFromGFF("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Counts/Ensembl/Homo_sapiens.GRCh37.75.gtf",
	# format="gtf")
# saveDb(TranscriptDb, file="/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/ensembl_v75_txdb_update.sqlite")


# DER for figure, ZSCAN23
pdf("test_SUGP1_structure.pdf",h=3.5,w=7)
gr = GRanges("chr19", IRanges(19387300,19431350))
# par(mar=c(1,5,1,1))
# layout(matrix(c(1,1, 2,2,2), nc= 1))
# meanCov = import("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/Adult_coverage.bw",
	# which = gr)
# plot(log2(meanCov$score+1) ~ start(meanCov),type="l",xaxt="n", xlab="",
	# ylab="Mean Cov",cex.axis=1.4,cex.lab=1.5)
# abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$SNP=="rs1233578"], ":",2)[1],"-",1)),
	# as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$SNP=="rs1233578"],"-",2)[1], "\\("))),lty=2)
par(mar=c(5,5,1,1))
plotTranscripts(gr, txByExon)
# abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$SNP=="rs1233578"], ":",2)[1],"-",1)),
	# as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$SNP=="rs1233578"],"-",2)[1], "\\("))),lty=2)
dev.off()


# # jxn for figure
# # chr5:137,837,170-138,284,761
# gr2 = GRanges("chr5", IRanges(137837170,138280000))

# pdf("plots/CTNNA1_structure.pdf",h=3.5,w=7)
# par(mar=c(5,5,1,1))
# plotTranscripts(gr2, txByExon)
# rect(xleft = 137946779, xright = 137988752,
	# ytop = 1, ybottom = -1,col="red")
# dev.off()

# # der for figure
# # chr3:180,730,454-181,525,311
# gr = GRanges("chr3", IRanges(180730454,181525311))
# meanCov = import("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/Adult_coverage.bw",
	# which = gr)

# pdf("plots/SOX2OT_structure.pdf",h=3.5,w=7)
# par(mar=c(1,5,1,1))
# layout(matrix(c(1,1, 2,2,2), nc= 1))
# plot(log2(meanCov$score+1) ~ end(meanCov),type="l",xaxt="n", xlab="",
	# ylab="Mean Cov",cex.axis=1.4,cex.lab=1.5, xlim = c(start(gr), end(gr)))
# abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"], ":",2)[1],"-",1)),
	# as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"],"-",2)[1], "\\("))),lty=2)
# par(mar=c(5,5,1,1))
# plotTranscripts(gr, txByExon)
# abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"], ":",2)[1],"-",1)),
	# as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"],"-",2)[1], "\\("))),lty=2)
# dev.off()


# #######

	
# ## check zhu
# gg = c("SF3B1", "PCCB", "ENDOG", "ZDHHC5", "SNX19", "C12orf24",
	# "ABCB9", "BAG5", "PSMA4","NMB", "NMRAL1", "KCTD13", "DUS2L",
	# "C17ORF39", "GATAD2A" ,"IRF3")
# gg[which(gg%in% pgcEqtlsFinal$Symbol)]
# gg[which( ! gg%in% pgcEqtlsFinal$Symbol)]


# ## devel stats
# load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
# load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")
# statListPgc = endoapply(statList, function(x) x[names(x) %in% pgcEqtlsFinal$Feature])

# devStatsPgc = unlist(statListPgc)
# names(devStatsPgc) = ss(names(devStatsPgc),"\\.",2)
# devStatsPgc = devStatsPgc[pgcEqtlsFinal$Feature]

# ## tabulate
# mean(devStatsPgc$p_bonf < 0.05, na.rm=TRUE) # 58.7%
# sum(sapply(statList, function(x) sum(x$p_bonf < 0.05,na.rm=TRUE)))/sum(lengths(statList)) # 16%
# devTable = table(pgcEqtlsFinal$riskDir[devStatsPgc$p_bonf < 0.05], 
	# ifelse(devStatsPgc$ageCorr[devStatsPgc$p_bonf < 0.05] < 0, "Fetal","Postnatal"))
# chisq.test(devTable)
# getOR(devTable)

# #####
# ## check isoswitch
# isoGenes = unique(unlist(lapply(switchList, rownames)))
# isoGenes = isoGenes[!is.na(isoGenes) & !grepl("-", isoGenes)]
# table(unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)

# ## whats bg?
# bgEns = unique(unlist(lapply(statList,
	# function(x) x$EnsemblGeneID[which(x$p_bonf < 0.05)])))
# bgEns = bgEns[!is.na(bgEns) & !grepl("-", bgEns)]

# # any
# topleft = sum(unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)
# topright = sum(! unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)
# bottomleft = length(unique(isoGenes)) - topleft
# bottomright = length(bgEns) - topleft - topright - bottomleft
# isoTab= matrix(c(topleft, topright, bottomleft, bottomright),
	# nrow = 2, byrow = TRUE,dimnames=list(c("pgcEqtl","noPgcEql"), 
		# c("isoSwitch", "noSwitch")))
# chisq.test(isoTab)$p.value
# getOR(isoTab)

# # cond indep
# topleft = sum(unique(pgcEqtlsIndept$EnsemblID) %in% isoGenes)
# topright = sum(! unique(pgcEqtlsIndept$EnsemblID) %in% isoGenes)
# bottomleft = length(unique(isoGenes)) - topleft
# bottomright = length(bgEns) - topleft - topright - bottomleft
# isoTab_cond= matrix(c(topleft, topright, bottomleft, bottomright),
	# nrow = 2, byrow = TRUE,dimnames=list(c("pgcEqtl","noPgcEql"), 
		# c("isoSwitch", "noSwitch")))
# chisq.test(isoTab_cond)$p.value
# getOR(isoTab_cond)
# prop.table(isoTab_cond,1)
