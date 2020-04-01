#
library('SummarizedExperiment')
library('getopt')
library('limma')
library('edgeR')
library('devtools')
library('jaffelab')
library('lmerTest')

## SNP / features to test (from PZ)
amyg = read.csv("raggr_881_leadvarsigfeat_amyg.csv", row.names=1)
amyg$Region = "Amygdala"
sacc = read.csv("raggr_881_leadvarsigfeat_sacc.csv", row.names=1)
sacc$Region = "sACC"
dlpfc = read.csv("raggr_881_leadvarsigfeat_dlpfc.csv", row.names=1)
dlpfc$Region = "DLPFC"
colnames(dlpfc)[6] = "Type"
# combine to one df
leadVarTests = rbind(amyg,sacc,dlpfc)
leadVarTests = leadVarTests[,c("Region","SNP","gene","Type","Symbol","leadVariant","hg19POS")]
colnames(leadVarTests)[3] = "feature"
leadVarTests$feature = as.character(leadVarTests$feature)



## combined data (pca, pd, snps)
load("3regions_expressionPCs.rda", verbose=TRUE)
load("3regions_combined_data.rda", verbose=TRUE)
pd = pdCombined
pd$BrainRegion = relevel(pd$BrainRegion, ref="sACC")

## fix jxn names
leadVarTests$feature = gsub("\\(\\-","\\(\\*", leadVarTests$feature)
leadVarTests$feature = gsub("\\(\\+","\\(\\*", leadVarTests$feature)

## fix exon names
leadVarTests$feature2 = NA
leadVarTests$featureE = NA
leadVarTests$feature2 = exonMapD$exonPos[match(leadVarTests$feature, exonMapD$exon_libdID)]  ## with e#### ids
leadVarTests$featureE = exonMapD$exonPos[match(leadVarTests$feature, exonMapD$exon_gencodeID)]  ## with ENSE### ids
# combine into feature2 (based on coord)
leadVarTests$feature2[is.na(leadVarTests$feature2)] = leadVarTests$feature[is.na(leadVarTests$feature2)]
leadVarTests$feature2[!is.na(leadVarTests$featureE)] = leadVarTests$featureE[!is.na(leadVarTests$featureE)] 

### Drop missing results (novel junctions not present in DLPFC)
dropInd = which(! leadVarTests$feature2 %in% rownames(rpkmCombined))
# leadVarTests$feature2[dropInd]
# [1] "chr4:94875901-94985640(*)"   "chr4:122190080-122196634(*)" "chr16:69336713-69339024(*)"  
# "chr3:167633533-167635696(*)" "chr4:94875901-94985640(*)"   "chr4:118111680-118114805(*)"
# "chr4:122190080-122196634(*)" "chr16:69336713-69339024(*)"

leadVarTests = leadVarTests[-dropInd,]

## drop duplicated rows (ones pz listed in multiple regions)
pair = paste0(leadVarTests$hg19POS, "_", leadVarTests$feature2)
dupInd = which(duplicated(pair))
leadVarTests = leadVarTests[-dupInd,]


## model
# Y = a + b1snp + b2brain(indicator_amyg) + b3brain(indicator_dlpfc) + b4snp*brain(amyg-vs-sacc) + b5snp*(dlpfc-vs-sacc) + covariates

n = nrow(leadVarTests)

sumList = list()
interDF = data.frame(feat=rep(NA,n), Symbol=rep(NA,n), leadVar=rep(NA,n), sampleSize=rep(NA,n), snpXamyg_pval=rep(NA,n), snpXDLPFC_pval=rep(NA,n))

pdf("example_interactions.pdf",w=13)
par(mar = c(5,6,2,2), cex.axis=1.1,cex.lab=1.8)
for(i in 1:n) {
	if(i %% 10 == 0) cat(".")
	# lv_i = leadVarTests$hg19POS[i]     ## identify by coords.
	lv_i = as.character(leadVarTests$hg19POS[i])		## previous line uses wrong index bc of factor
	snp_i = as.numeric(snpCombined[lv_i,]) 
	feat_i = leadVarTests$feature2[i]
	
	if (sum(is.na(snp_i))>0) {  # remove NA samples
		dropInd = which(is.na(snp_i))
		pdtemp = pd[-dropInd,]
		snptmp = snp_i[-dropInd]
		mdstmp = as.matrix(mdsCombined[-dropInd,1:5])
		pcatmp = pcaExp$x[-dropInd,1:5]
		rpkmtmp = rpkmCombined[,-dropInd]
	} else {
		pdtemp = pd
		snptmp = snp_i
		mdstmp = as.matrix(mdsCombined[,1:5])
		pcatmp = pcaExp$x[,1:5]
		rpkmtmp = rpkmCombined
	}
	sampsize = length(snptmp)

	# mod_i = model.matrix(~0 + snptmp + BrainRegion + snptmp*BrainRegion + AgeDeath + Sex + mdstmp + pcatmp, data=pdtemp)
	# mod_i = mod_i[,-2] ## remove sACC column
	# fit = lm(rpkmtmp[feat_i,] ~ mod_i)
	
	## aej edits, adding log2 xform
	fit = lmer(log2(rpkmtmp[feat_i,]+1) ~ snptmp + pdtemp$BrainRegion + 
		snptmp*pdtemp$BrainRegion + pdtemp$AgeDeath + pdtemp$Sex + mdstmp + pcatmp +
		(1|pdtemp$BrNum))
		
	sumList[[i]] = summary(fit)
	interDF$feat[i] = feat_i
	interDF$Symbol[i] = as.character(leadVarTests$Symbol[i])
	interDF$leadVar[i] = as.character(leadVarTests$SNP[i])
	interDF$sampleSize[i] = length(summary(fit)$res)
	interDF$snpXamyg_pval[i] = summary(fit)$coefficients[17,5]
	interDF$snpXDLPFC_pval[i] = summary(fit)$coefficients[18,5]
	interDF$snpRegionAnova[i] = as.data.frame(anova(fit))[7,6]
	
	## add plot
	boxplot(log2(rpkmtmp[feat_i,]+1) ~ snptmp*pdtemp$BrainRegion,
		xlab = as.character(leadVarTests$SNP[i]),
		ylab = paste0("log2 ", feat_i))
	legend("topleft", paste0("intANOVA p=",
		signif(interDF$snpRegionAnova[i],3)), cex=1.5)
}

dev.off()
# write.csv(interDF, file="bipolar_leadvar_interactions.csv")
# write.csv(interDF, file="bipolar_leadvar_interactions_lmer_inclAnova.csv")
write.csv(interDF, file="bipolar_leadvar_interactions_lmer_inclAnova_char.csv")

# interDF$qval_amyg = p.adjust(interDF$snpXamyg_pval, "fdr")
# interDF$qval_dlpfc = p.adjust(interDF$snpXDLPFC_pval, "fdr")

# interDF$bonf_amyg = p.adjust(interDF$snpXamyg_pval, "bonferroni")
# interDF$bonf_dlpfc = p.adjust(interDF$snpXDLPFC_pval, "bonferroni")

#####################
#### linear flavor

n = nrow(leadVarTests)

sumListLm = list()
interDFlm = data.frame(feat=rep(NA,n), Symbol=rep(NA,n), leadVar=rep(NA,n), sampleSize=rep(NA,n), snpXamyg_pval=rep(NA,n), snpXDLPFC_pval=rep(NA,n))

for(i in 1:n) {
	if(i %% 10 == 0) cat(".")
	# lv_i = leadVarTests$hg19POS[i]     ## identify by coords.
	lv_i = as.character(leadVarTests$hg19POS[i])		## previous line uses wrong index bc of factor
	snp_i = as.numeric(snpCombined[lv_i,]) 
	feat_i = leadVarTests$feature2[i]
	
	pdtemp = pd
	snptmp = snp_i
	mdstmp = as.matrix(mdsCombined[,1:5])
	pcatmp = pcaExp$x[,1:5]
	rpkmtmp = rpkmCombined

	## aej edits
	fit = lm(rpkmtmp[feat_i,] ~ snptmp + pdtemp$BrainRegion + 
		snptmp*pdtemp$BrainRegion + pdtemp$AgeDeath + pdtemp$Sex + mdstmp + pcatmp)
		
	sumListLm[[i]] = summary(fit)
	interDFlm$feat[i] = feat_i
	interDFlm$Symbol[i] = as.character(leadVarTests$Symbol[i])
	interDFlm$leadVar[i] = as.character(leadVarTests$SNP[i])
	interDFlm$sampleSize[i] = length(fit$fitted)
	interDFlm$snpXamyg_pval[i] = summary(fit)$coefficients[17,4]
	interDFlm$snpXDLPFC_pval[i] = summary(fit)$coefficients[18,4]
	interDFlm$snpRegionAnova[i] = as.data.frame(anova(fit))[7,5]
	
}

write.csv(interDFlm, file="bipolar_leadvar_interactions_linear_inclAnova_char.csv")

