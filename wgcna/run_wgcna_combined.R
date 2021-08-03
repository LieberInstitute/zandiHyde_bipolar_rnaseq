####### 
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)

## multithread
allowWGCNAThreads(8)

dir.create("rdas")

## load data
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/degradation_rse_BipSeq_BothRegions.rda")

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
				
## add ancestry 
load("../genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

###########
# filter ##
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

##########
## model #
##########
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 18
qSV_mat = prcomp(t(degExprs))$x[,1:k]

## join and move around region, dx and interaction for cleaning
modQsva = cbind(modJoint[,c(1:4,14,5:13)], qSV_mat)

## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)

#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,
                               networkType = "signed", verbose = 5)
cat(sftthresh1$powerEstimate)
save(sftthresh1, file = "rdas/power_object.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "rdas/wgcna_signed_TOM")
fNames = rownames(geneExprs)
save(net, fNames, file = "rdas/constructed_network_signed_bicor.rda")

########################
## by region
modRegion =  model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))
rIndexes = splitit(rse_gene$BrainRegion)

## clean by region
geneExprs_list = mclapply(rIndexes, function(ii) {
	degExprs = log2(assays(cov_rse[,ii])$count+1)
	k = num.sv(degExprs, modRegion[ii,]) 
	m = cbind(modRegion[ii,], prcomp(t(degExprs))$x[,1:k])
	cleaningY(geneExprs[,ii], m, P=3)
},mc.cores=2)

## threshold
thresh_list = mclapply(geneExprs_list, function(y) {
	pickSoftThreshold(t(y), powerVector = powers,
                               networkType = "signed", verbose = 5)
},mc.cores=2)

## networks
net_list = lapply(1:2, function(i) {
	blockwiseModules(t(geneExprs_list[[i]]), power = thresh_list[[i]]$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = paste0("rdas/wgcna_signed_TOM_region",
								names(rIndexes)[i]))
})
save(net_list, net, fNames, file = "rdas/constructed_network_signed_bicor.rda")

