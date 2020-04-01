#### 
library(GenomicFeatures)
library('BSgenome.Hsapiens.UCSC.hg38')
library(derfinderPlot)
library(jaffelab)
library(SummarizedExperiment)
library(RColorBrewer)


# ########################################################
# ############## tx plots

# ## GENCODE
# chrInfo = getChromInfoFromUCSC("hg38")
# chrInfo$chrom = as.character(chrInfo$chrom)
# chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
# chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
# si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

# gencode_v25 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
# seqlevels(gencode_v25, pruning.mode="coarse") = paste0("chr", c(1:22,"X","Y"))
# seqinfo(gencode_v25) = si
# gencode_v25_txdb = makeTxDbFromGRanges(gencode_v25)
# txByExon = exonsBy(gencode_v25_txdb, use.names=TRUE)
# geneMap = as.data.frame(gencode_v25[which(gencode_v25$type == "gene"),])

# amyg = read.csv("../conditional/raggr_881_snps_amyg_eqtls_fdr01.csv", row.names=1)
# dlp = read.csv("../conditional/raggr_881_snps_dlpfc_eqtls_fdr01.csv")
# sacc = read.csv("../conditional/raggr_881_snps_sacc_eqtls_fdr01.csv")

# ## Add leadVariant position
# load("../../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
# snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

# riskLoci = read.csv("../rAggr_results_881.csv", stringsAsFactors=FALSE)	# 13,592 snps
# colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
# riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
# riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

# riskLoci$hg38POS1 = snpMap$pos_hg38[match(riskLoci$hg19POS1, snpMap$pos_hg19)]
# riskLoci$hg38POS2 = snpMap$pos_hg38[match(riskLoci$hg19POS2, snpMap$pos_hg19)]

# ## a couple manually, from dbSNP
# riskLoci$hg38POS1[which(riskLoci$SNP1_Name == "rs10035291:80796368:T:C")] = 81500549
# riskLoci$hg38POS2[which(riskLoci$SNP2_Name == "rs10035291:80796368:T:C")] = 81500549
# riskLoci$hg38POS1[which(riskLoci$SNP1_Name == "rs112114764:42201041:G:T")] = 44123673
# riskLoci$hg38POS2[which(riskLoci$SNP2_Name == "rs112114764:42201041:G:T")] = 44123673
# riskLoci$hg38POS1[which(riskLoci$SNP1_Name == "rs59134449:111745562:G:GT")] = 109985804
# riskLoci$hg38POS2[which(riskLoci$SNP2_Name == "rs59134449:111745562:G:GT")] = 109985804


# amyg$leadVariant_hg19POS = riskLoci$hg19POS2[match(amyg$leadVariant, riskLoci$SNP2_Name)]
# amyg$leadVariant_hg38POS = riskLoci$hg38POS2[match(amyg$leadVariant, riskLoci$SNP2_Name)]
# amyg$IndexSNP_hg38POS = riskLoci$hg38POS1[match(amyg$IndexSNP_hg19POS, riskLoci$hg19POS1)]

# dlp$leadVariant_hg19POS = riskLoci$hg19POS2[match(dlp$leadVariant, riskLoci$SNP2_Name)]
# dlp$leadVariant_hg38POS = riskLoci$hg38POS2[match(dlp$leadVariant, riskLoci$SNP2_Name)]
# dlp$IndexSNP_hg38POS = riskLoci$hg38POS1[match(dlp$IndexSNP_hg19POS, riskLoci$hg19POS1)]

# sacc$leadVariant_hg19POS = riskLoci$hg19POS2[match(sacc$leadVariant, riskLoci$SNP2_Name)]
# sacc$leadVariant_hg38POS = riskLoci$hg38POS2[match(sacc$leadVariant, riskLoci$SNP2_Name)]
# sacc$IndexSNP_hg38POS = riskLoci$hg38POS1[match(sacc$IndexSNP_hg19POS, riskLoci$hg19POS1)]


# ## genes to plot:
# gen = c("SSBP2","PACS1","SUGP1","GATAD2A","ASB16-AS1","ZCCHC2","GRIN2A","CSPG4P12","SCN2A","NEK4","GNL3",
		# "LRRC57","GANC","LMAN2L","ADD3","STK4","CTSF","HIST2H4A","TRANK1","LRRFIP2")
# indsnp = c("rs10035291:80796368:T:C","rs10896090:65945186:A:G","rs111444407:19358207:C:T","rs111444407:19358207:C:T",
		# "rs112114764:42201041:G:T","rs11557713:60243876:G:A","rs11647445:9926966:T:G","rs139221256:85357857:T:TA",
		# "rs17183814:166152389:G:A","rs2302417:52814256:T:A","rs2302417:52814256:T:A","rs4447398:42904904:A:C",
		# "rs4447398:42904904:A:C","rs57195239:97376407:A:AT","rs59134449:111745562:G:GT","rs6130764:43750410:T:C",
		# "rs7122539:66662731:G:A","rs7544145:150138699:C:T","rs9834970:36856030:T:C","rs9834970:36856030:T:C")
# leadsnp = c("rs6887937:80748607:A:G","rs6503490:42246656:A:G","rs7080576:111710817:T:C")


		
# geneInfo = data.frame(symbol=gen, 
			# chr=geneMap$seqnames[match(gen, geneMap$gene_name)],
			# start=geneMap$start[match(gen, geneMap$gene_name)],
			# end=geneMap$end[match(gen, geneMap$gene_name)],
			# strand=geneMap$strand[match(gen, geneMap$gene_name)],
			# indexsnp=indsnp)
# geneInfo$indexsnp_pos = riskLoci$hg38POS2[match(indsnp, riskLoci$SNP2_Name)]

# leadsnp=NA
# geneInfo$leadsnp_pos = NA
# geneInfo$leadsnp[which(geneInfo$symbol %in% c("SSBP2","ASB16-AS1","ADD3"))] = leadsnp
# geneInfo$leadsnp_pos[which(geneInfo$symbol %in% c("SSBP2","ASB16-AS1","ADD3"))] = riskLoci$hg38POS2[match(leadsnp, riskLoci$SNP2_Name)]
			
# pdf("eqtl_plots/tx_structure.pdf",h=3.5,w=10)
# par(mar=c(5,3,2,2))

# for (i in 1:nrow(geneInfo)) {
	# ## select gene coords
	# chrom = geneInfo$chr[i]
	# minRange = geneInfo$start[i]
	# maxRange = geneInfo$end[i]
	# strnd=geneInfo$strand[i]
	# sym=as.character(geneInfo$symbol[i])
	
	# if (sym == "PACS1") { 
		# maxRange=maxRange+1500
		# strnd="*"
		# sym="PACS1 / RP11-755F10.1"
	# }

	# gr = GRanges(chrom, IRanges(minRange,maxRange), strnd)
	# ov <- findOverlaps(gr, txByExon)
	# txList <- split(txByExon[subjectHits(ov)], queryHits(ov))
	# poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
	# yrange <- range(poly.data$exon$y)
	# xrange <- c(start(gr),end(gr))

	# plot(0,0, type="n", xlim = xrange , ylim = yrange + c(-0.75, 0.75),
		# yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
		# cex.axis = 1.3, cex.lab =1.3, cex.main =1.3, main=sym)
	# # rect(xleft=66220288, xright=66220791, ybottom=-100, ytop=100, border=FALSE, col="peachpuff")
	# # box(which = "plot") ## fix border
	
	# # if (!is.na(geneInfo$leadsnp[i])) { ## leadVar is diff than index
	# # abline(v=geneInfo$indexsnp_pos[i], lty=2, col="ivory4")  	## index
	# # abline(v=geneInfo$leadsnp_pos[i], lty=2, col="maroon")  	## leadVar
	# # } else { ## leadVar and index are same
	# # abline(v=geneInfo$indexsnp_pos[i], lty=2, col="maroon")  	## index
	# # }
	
	# derfinderPlot:::.plotPolygon(poly.data$exon, 'blue')
	# derfinderPlot:::.plotPolygon(poly.data$intron, 'lightblue')
	# yTx <- unique(yrange / abs(yrange))
	# if(length(yTx) > 1) {
		# axis(2, yTx, c('-', '+')[c(-1, 1) %in% yTx], 
			# tick = FALSE, las = 1, cex.axis = 3)
		# abline(h = 0, lty = 3)
	# }

# }
# dev.off()
	

	
	
	
########################################################
############## eQTL plots

gen = c("SSBP2","PACS1","SUGP1","GATAD2A","ASB16-AS1","ZCCHC2","GRIN2A","CSPG4P12","SCN2A","NEK4","GNL3",
		"LRRC57","GANC","LMAN2L","ADD3","STK4","CTSF","HIST2H4A","TRANK1","LRRFIP2","XPNPEP1")
		
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_amygdala.rda")
geneMap = as.data.frame(rowRanges(rse_gene))
geneMap = geneMap[geneMap$Symbol %in% gen , c("gencodeID","Symbol")]
exonMap = as.data.frame(rowRanges(rse_exon))
exonMap = exonMap[exonMap$Symbol %in% gen , c("gencodeID","Symbol")]
jxnMap = as.data.frame(rowRanges(rse_jxn))
jxnMap = jxnMap[jxnMap$newGeneSymbol %in% gen , c("newGeneID","newGeneSymbol")]
names(jxnMap) = c("gencodeID","Symbol")
txMap = as.data.frame(rowRanges(rse_tx))
txMap = txMap[txMap$gene_name %in% gen , c("gene_id","gene_name")]
names(txMap) = c("gencodeID","Symbol")
map_amyg = rbind(geneMap, exonMap, jxnMap, txMap)
map_amyg$Type = c(rep("Gene",nrow(geneMap)),rep("Exon",nrow(exonMap)),rep("Jxn",nrow(jxnMap)),rep("Tx",nrow(txMap)))

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_sacc.rda")
geneMap = as.data.frame(rowRanges(rse_gene))
geneMap = geneMap[geneMap$Symbol %in% gen , c("gencodeID","Symbol")]
exonMap = as.data.frame(rowRanges(rse_exon))
exonMap = exonMap[exonMap$Symbol %in% gen , c("gencodeID","Symbol")]
jxnMap = as.data.frame(rowRanges(rse_jxn))
jxnMap = jxnMap[jxnMap$newGeneSymbol %in% gen , c("newGeneID","newGeneSymbol")]
names(jxnMap) = c("gencodeID","Symbol")
txMap = as.data.frame(rowRanges(rse_tx))
txMap = txMap[txMap$gene_name %in% gen , c("gene_id","gene_name")]
names(txMap) = c("gencodeID","Symbol")
map_sacc = rbind(geneMap, exonMap, jxnMap, txMap)
map_sacc$Type = c(rep("Gene",nrow(geneMap)),rep("Exon",nrow(exonMap)),rep("Jxn",nrow(jxnMap)),rep("Tx",nrow(txMap)))

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
geneMap = as.data.frame(rowRanges(rse_gene))
geneMap = geneMap[geneMap$Symbol %in% gen , c("gencodeID","Symbol")]
exonMap = as.data.frame(rowRanges(rse_exon))
exonMapIDs = exonMap[,c("exon_gencodeID","exon_libdID")]
rownames(exonMap) = exonMap$exon_libdID
exonMap = exonMap[exonMap$Symbol %in% gen , c("gencodeID","Symbol")]
jxnMap = as.data.frame(rowRanges(rse_jxn))
jxnMap = jxnMap[jxnMap$newGeneSymbol %in% gen , c("newGeneID","newGeneSymbol")]
names(jxnMap) = c("gencodeID","Symbol")
txMap = as.data.frame(rowRanges(rse_tx))
txMap = txMap[txMap$gene_name %in% gen , c("gene_id","gene_name")]
names(txMap) = c("gencodeID","Symbol")
map_dlpfc = rbind(geneMap, exonMap, jxnMap, txMap)
map_dlpfc$Type = c(rep("Gene",nrow(geneMap)),rep("Exon",nrow(exonMap)),rep("Jxn",nrow(jxnMap)),rep("Tx",nrow(txMap)))

## make jxns match exprs (since dlpfc is *)
rownames(map_amyg) = ss(rownames(map_amyg), "\\(")
rownames(map_sacc) = ss(rownames(map_sacc), "\\(")
rownames(map_dlpfc) = ss(rownames(map_dlpfc), "\\(")

## load expression
load("residualized_exprs_3regions.rda", verbose=TRUE)
rownames(amyg_log2exprs_clean) = ss(rownames(amyg_log2exprs_clean), "\\(")
rownames(sacc_log2exprs_clean) = ss(rownames(sacc_log2exprs_clean), "\\(")
rownames(dlpfc_log2exprs_clean) = ss(rownames(dlpfc_log2exprs_clean), "\\(")

## change exon row names
stopifnot(identical(rownames(dlpfc_log2exprs_clean)[25137:421954], exonMapIDs$exon_gencodeID))
rownames(dlpfc_log2exprs_clean)[25137:421954] = exonMapIDs$exon_libdID


## load genotypes
load("../../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda", verbose=TRUE)
snpMap_amyg = snpMap
snpMap_amyg$snp2 = ss(snpMap_amyg$SNP,":")
snp_amyg = snp[,pd_amyg$BrNum]
colnames(snp_amyg) = pd_amyg$RNum
rownames(snp_amyg) = ss(rownames(snp_amyg),":")

snp_sacc = snp[,pd_sacc$BrNum]
colnames(snp_sacc) = pd_sacc$RNum
rownames(snp_sacc) = ss(rownames(snp_sacc),":")

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda", verbose=TRUE)
snpMap_dlpfc = snpMap
snpMap_dlpfc$snp2 = ss(snpMap_dlpfc$SNP,":")
snp_dlpfc = snp[,pd_dlpfc$BrNum]
colnames(snp_dlpfc) = pd_dlpfc$RNum
rownames(snp_dlpfc) = ss(rownames(snp_dlpfc),":")

### plot expression by snp for combos given by Peter
# exprs by indexsnp - amyg,sacc,dlpfc
# exprs by leadsnp - amyg,sacc,dlpfc

toplot = read.csv("pz_toplot.csv", stringsAsFactors=FALSE)
toplot$indata = toplot$snpi2 %in% rownames(snp_amyg)  ## only missing index snps are missing
toplot = toplot[toplot$indata,]

snpMap_toplot_amyg = snpMap_amyg[match(toplot$snpi2,snpMap_amyg$snp2),]
snpMap_toplot_dlpfc = snpMap_dlpfc[match(toplot$snpi2,snpMap_dlpfc$snp2),]
identical(snpMap_toplot_amyg, snpMap_toplot_dlpfc)  ## TRUE

toplot$ref = snpMap_toplot_amyg$newRef
toplot$count = snpMap_toplot_amyg$newCount

## PLOTS
for (i in 1:nrow(toplot)) {
# for (i in c(1,2,4,6,7,9,10,12:24)) {
symi = toplot$symi[i]
snpi = toplot$snpi2[i]
refi = toplot$ref[i]
counti = toplot$count[i]
xaxisi = c(paste0(refi,refi),paste0(refi,counti),paste0(counti,counti))
print(paste0(symi, ", ",snpi))

## go through all features of that sym

pdf(paste0("eqtl_plots_labeled/4feat_eqtl_3region_",symi,"_",snpi,"_2.pdf"), h=4, w=9)
par(mfrow=c(1,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))	

allfeat1 = map_amyg[which(map_amyg$Symbol == symi),]
allfeat2 = map_sacc[which(map_sacc$Symbol == symi),]
allfeat3 = map_dlpfc[which(map_dlpfc$Symbol == symi),]
## only plot feats present in all 3 (drop novel jxns only present in 1 or 2)
allfeat = allfeat1[which(rownames(allfeat1) %in% rownames(allfeat3)),]

for (f in 1:nrow(allfeat)) {
	feati = rownames(allfeat)[f]
	typei = allfeat$Type[f]		
	## plot
		boxplot(amyg_log2exprs_clean[feati,] ~ snp_amyg[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("Amygdala\n",symi,"\n",feati," (",typei,")"), 
				ylim = c(range(amyg_log2exprs_clean[feati,])), outline=FALSE,
				xaxt="n")
		axis(1, at=1:3, labels=xaxisi)
		points(amyg_log2exprs_clean[feati,] ~ jitter(snp_amyg[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_amyg[snpi,])+2,cex=1.5)
		# legend("top",paste0("p=",p_i))
		
		boxplot(sacc_log2exprs_clean[feati,] ~ snp_sacc[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("sACC\n",symi,"\n",feati," (",typei,")"), 
				ylim = c(range(sacc_log2exprs_clean[feati,])), outline=FALSE,
				xaxt="n")
		axis(1, at=1:3, labels=xaxisi)
		points(sacc_log2exprs_clean[feati,] ~ jitter(snp_sacc[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_sacc[snpi,])+2,cex=1.5)
				   
		boxplot(dlpfc_log2exprs_clean[feati,] ~ snp_dlpfc[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("DLPFC\n",symi,"\n",feati," (",typei,")"), 
				ylim = c(range(dlpfc_log2exprs_clean[feati,])), outline=FALSE,
				xaxt="n")
		axis(1, at=1:3, labels=xaxisi)
		points(dlpfc_log2exprs_clean[feati,] ~ jitter(snp_dlpfc[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_dlpfc[snpi,])+2,cex=1.5)
	}
	dev.off()
}








### SCN2A individual plot:

i = 10	# Only SCN2A
symi = toplot$symi[i]
snpi = toplot$snpi2[i]
refi = toplot$ref[i]
counti = toplot$count[i]
xaxisi = c(paste0(refi,refi),paste0(refi,counti),paste0(counti,counti))
print(paste0(symi, ", ",snpi))

## go through all features of that sym

pdf(paste0("eqtl_plots_consistentAxis/4feat_eqtl_3region_",symi,"_",snpi,"_2.pdf"), h=4, w=9)
par(mfrow=c(1,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))	

allfeat1 = map_amyg[which(map_amyg$Symbol == symi),]
allfeat2 = map_sacc[which(map_sacc$Symbol == symi),]
allfeat3 = map_dlpfc[which(map_dlpfc$Symbol == symi),]
## only plot feats present in all 3 (drop novel jxns only present in 1 or 2)
allfeat = allfeat1[which(rownames(allfeat1) %in% rownames(allfeat3)),]

f = 70	# Only the one jxn
	feati = rownames(allfeat)[f]
	typei = allfeat$Type[f]		
	## plot
	snpline_amyg = snp_amyg[snpi,]+1
		boxplot(amyg_log2exprs_clean[feati,] ~ snp_amyg[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("Amygdala\n",symi,"\n",feati," (",typei,")"), 
				ylim = c(-.5,5), outline=FALSE,
				xaxt="n")
		axis(1, at=1:3, labels=xaxisi)
		abline(lm(amyg_log2exprs_clean[feati,] ~ snpline_amyg), col="lightblue", lwd=1.25)
		points(amyg_log2exprs_clean[feati,] ~ jitter(snp_amyg[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_amyg[snpi,])+2,cex=1.5)
		# legend("top",paste0("p=",p_i))	
	snpline_sacc = snp_sacc[snpi,]+1
	snpline_dlp = snp_dlpfc[snpi,]+1
	exprsline_sacc = sacc_log2exprs_clean[feati,]-4.4		# to continue best fit into other plot
	exprsline_dlp = dlpfc_log2exprs_clean[feati,]+2.85
		boxplot(sacc_log2exprs_clean[feati,] ~ snp_sacc[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("sACC\n",symi,"\n",feati," (",typei,")"), 
				ylim = c(-.5,5), outline=FALSE,
				xaxt="n")
		axis(1, at=1:3, labels=xaxisi)
		abline(lm(sacc_log2exprs_clean[feati,] ~ snpline_sacc), col="lightblue", lwd=1.25)
		abline(lm(exprsline_dlp ~ snpline_dlp), col="lightgreen", lwd=1.25)
		points(sacc_log2exprs_clean[feati,] ~ jitter(snp_sacc[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_sacc[snpi,])+2,cex=1.5)
				   
		boxplot(dlpfc_log2exprs_clean[feati,] ~ snp_dlpfc[snpi,],
				xlab=snpi, ylab="Residualized Expression", 
				main=paste0("DLPFC\n",symi,"\n",feati," (",typei,")"), 
				xlim = c(0.5,3.5),
				ylim = c(-.5,5), outline=FALSE,
				xaxt="n")
		abline(lm(exprsline_sacc ~ snpline_sacc), col="lightblue", lwd=1.25)
		abline(lm(dlpfc_log2exprs_clean[feati,] ~ snpline_dlp), col="lightgreen", lwd=1.25)
		axis(1, at=1:3, labels=xaxisi)
		points(dlpfc_log2exprs_clean[feati,] ~ jitter(snp_dlpfc[snpi,]+1),
				   pch=21, 
				   bg=as.numeric(snp_dlpfc[snpi,])+2,cex=1.5)
	}
	dev.off()
	
	
	
	
	
	
	
	
	
	