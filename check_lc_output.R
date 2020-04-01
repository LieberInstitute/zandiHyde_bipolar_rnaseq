##
library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(sva)
## read in leafcutter output, case-control
fn = list.files("leafcutter_output", pattern = "_ds_", full=TRUE)
datList = lapply(fn, read.delim,as.is=TRUE)
names(datList) = list.files("leafcutter_output",pattern = "_ds_")

## read in junction data

#### load in data
load("data/zandiHypde_bipolar_rseJxn_n511.rda")
load("data/zandiHypde_bipolar_rseGene_n511.rda")
load("data/degradation_rse_BipSeq_BothRegions.rda")

identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
				
## add ancestry 
load("genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## fix dx
rse_jxn$PrimaryDx[rse_jxn$PrimaryDx == "Other"] = "Bipolar"
rse_jxn$PrimaryDx = droplevels(rse_jxn$PrimaryDx)
rse_jxn$PrimaryDx = factor(rse_jxn$PrimaryDx, c("Control", "Bipolar"))

## and in gene level for modeling
rse_gene$PrimaryDx[rse_gene$PrimaryDx == "Other"] = "Bipolar"
rse_gene$PrimaryDx = droplevels(rse_gene$PrimaryDx)
rse_gene$PrimaryDx = factor(rse_gene$PrimaryDx, c("Control", "Bipolar"))

#### model
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))
## with qSVs
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 18
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 87%
modQsva = cbind(modJoint, qSV_mat)

##############################

## split up
amyg_intron = datList[[1]]
amyg_intron = amyg_intron[amyg_intron$intron != "intron",]
amyg_intron$deltapsi = as.numeric(amyg_intron$deltapsi)
amyg_intron$logef = as.numeric(amyg_intron$logef)
amyg_cluster = datList[[2]]
sacc_intron = datList[[3]]
sacc_intron = sacc_intron[sacc_intron$intron != "intron",]
sacc_intron$deltapsi = as.numeric(sacc_intron$deltapsi)
sacc_intron$logef = as.numeric(sacc_intron$logef)
sacc_cluster = datList[[4]]

# clena up
amyg_cluster$cluster_id = ss(amyg_cluster$cluster, ":", 2)
sacc_cluster$cluster_id = ss(sacc_cluster$cluster, ":", 2)
amyg_intron$cluster_id = ss(amyg_intron$intron, ":", 4)
sacc_intron$cluster_id = ss(sacc_intron$intron, ":", 4)

## match up
amyg_intron$rowMatch = match(amyg_intron$cluster_id, amyg_cluster$cluster_id)
sacc_intron$rowMatch = match(sacc_intron$cluster_id, sacc_cluster$cluster_id)

## add our junction ids
amyg_intron$jxnID = paste0(ss(amyg_intron$intron,":"), ":", 
	ss(amyg_intron$intron,":", 2), "-", 
	as.numeric(ss(amyg_intron$intron,":", 3))-1, "(",
	ss(amyg_intron$intron,"_", 3), ")")
sacc_intron$jxnID = paste0(ss(sacc_intron$intron,":"), ":", 
	ss(sacc_intron$intron,":", 2), "-", 
	as.numeric(ss(sacc_intron$intron,":", 3))-1, "(",
	ss(sacc_intron$intron,"_", 3), ")")
	
## match up there
amyg_intron$jxnMatch = match(amyg_intron$jxnID, rownames(rse_jxn))
sacc_intron$jxnMatch = match(sacc_intron$jxnID, rownames(rse_jxn))

## load results
load("case_control/bipolarControl_deStats_byRegion_qSVAjoint.rda")
amyg_stats = statOut[match(amyg_intron$jxnID,rownames(statOut)),]
sacc_stats = statOut[match(sacc_intron$jxnID,rownames(statOut)),]

## overall plots
plot(amyg_stats$logFC_Amyg, amyg_intron$deltapsi,pch=21, bg="grey")
abline(h=0,v=0,col="blue", lty=2)
plot(amyg_stats$logFC_Amyg, amyg_intron$logef,pch=21, bg="grey",
	ylim = c(-2,2), xlim = c(-2,2))
abline(h=0,v=0,col="blue", lty=2)
cor(amyg_stats$logFC_Amyg, amyg_intron$deltapsi)
abline(h=0,v=0,col="blue", lty=2)
cor(amyg_stats$logFC_Amyg, amyg_intron$logef)

plot(amyg_intron$deltapsi, amyg_intron$logef, ylim = c(-2,2))


##############################
### now lets spot check some #

rowData(rse_jxn)$Length = 100
jRp10m = recount::getRPKM(rse_jxn, "Length")
clean_exprs = 

###### cluster 1
amyg_cluster[1,]
x = amyg_intron[which(amyg_intron$rowMatch == 1),]
x
match(x$jxnID, sacc_intron$jxnID)

rowData(rse_jxn[x$jxnMatch,])

mypar(1,2)
y = jRp10m[x$jxnMatch,]
boxplot(log2(y[1,]+1) ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion,
	names = c("AMYG:CONT","AMYG:BPD", "sACC:CONT", "sACC:BPD"),
	ylim = c(0,14), xlab="", ylab = "log2(RP10M + 1)",
	main = rownames(y)[1])
boxplot(log2(y[2,]+1) ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion,
	names = c("AMYG:CONT","AMYG:BPD", "sACC:CONT", "sACC:BPD"),
	ylim = c(0,14), xlab="", ylab = "log2(RP10M + 1)",
	main = rownames(y)[2])

table(y[2,] == 0, paste0( rse_jxn$PrimaryDx , ":" , rse_jxn$BrainRegion))

### clean expression 
y_clean = cleaningY(log2(jRp10m[x$jxnMatch,]+1), modQsva, P=3)
boxplot(y_clean[1,] ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion,
	names = c("AMYG:CONT","AMYG:BPD", "sACC:CONT", "sACC:BPD"),
	ylim = c(0,12), xlab="", ylab = "Cleaned log2(RP10M + 1)",
	main = rownames(y)[1])
boxplot(y_clean[2,] ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion,
	names = c("AMYG:CONT","AMYG:BPD", "sACC:CONT", "sACC:BPD"),
	ylim = c(0,12), xlab="", ylab = "Cleaned log2(RP10M + 1)",
	main = rownames(y)[2])


## cluster 2
amyg_cluster[2,]
x = amyg_intron[which(amyg_intron$rowMatch == 2),]

match(x$jxnID, sacc_intron$jxnID)

rowData(rse_jxn[x$jxnMatch,])

y = jRp10m[x$jxnMatch,]
mypar(3,3)
for(i in 1:9) {
boxplot(log2(y[i,]+1) ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion)
}
boxplot(log2(y[2,]+1) ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion)
table(y[2,] == 0, paste0( rse_jxn$PrimaryDx , ":" , rse_jxn$BrainRegion))


## cluster 3
amyg_cluster[3,]
x = amyg_intron[which(amyg_intron$rowMatch == 3),]

match(x$jxnID, sacc_intron$jxnID)

rowData(rse_jxn[x$jxnMatch,])

y = jRp10m[x$jxnMatch,]
mypar(1,3)
for(i in 1:3) {
boxplot(log2(y[i,]+1) ~ rse_jxn$PrimaryDx * rse_jxn$BrainRegion)
}
table(y[2,] == 0, paste0( rse_jxn$PrimaryDx , ":" , rse_jxn$BrainRegion))

## top junctions?
jxnStats = statOut[order(statOut$P.Value_Amyg),]
jxnStats = jxnStats[grep("^chr", rownames(jxnStats)),]

amyg_intron[match(rownames(jxnStats)[1], amyg_intron$jxnID),]
amyg_cluster[amyg_cluster$cluster_id == "clu_10420_+",]
amyg_intron[amyg_intron$cluster_id == "clu_10420_+",]

########################
## check splice qQTLs

sqtls_sacc = read.delim("leafcutter_output/sacc.win.500000.bp.perm.0.cis.nominal.threaded.eqtl.allpairs.fdr.annotated.hg38gencodev25-fs20190926.txt",
	as.is=TRUE)
colnames(sqtls_sacc)[1:(ncol(sqtls_sacc)-1)] = colnames(sqtls_sacc)[2:(ncol(sqtls_sacc))]
sqtls_sacc = sqtls_sacc[,-ncol(sqtls_sacc)]
colnames(sqtls_sacc)[1] = "intron"

## add our junction ID
sqtls_sacc$jxnID = paste0("chr",ss(sqtls_sacc$intron,":"), ":", 
	ss(sqtls_sacc$intron,":", 2), "-", 
	as.numeric(ss(sqtls_sacc$intron,":", 3))-1, "(",
	ss(sqtls_sacc$intron,"_", 3), ")")
	
## make eQTL ID
sqtls_sacc$eqtlID = paste0(sqtls_sacc$snp, ";", sqtls_sacc$jxnID)

######################
## junction eQTLs
load("eqtl/genomewide/rdas/matrixEqtl_output_sacc_genomewide_jxn.rda")
eqtl_sacc = meJxn$cis$eqtl
eqtl_sacc$eqtlID = paste0(eqtl_sacc$snps, ";", eqtl_sacc$gene)

##################
## match up
mm = match(sqtls_sacc$eqtlID, eqtl_sacc$eqtlID)


## load SNP data
load("eqtl/genomewide/rdas/overlappingSNPs.rda")  # snpMapKeep
load("genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

snpInd = which(rownames(snpMap) %in% rownames(snpMapKeep) & !is.na(snpMap$pos_hg38))
snpMap = snpMap[snpInd,]
snp = snp[snpInd,]

### look at unique hits
sqtls_sacc_uniq = sqtls_sacc[!duplicated(sqtls_sacc$variant_id),]
