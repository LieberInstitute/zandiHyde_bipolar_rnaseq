##
library(jaffelab)
library(SummarizedExperiment)
library(LIBDpheno)

# phenotype data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda")
pd = colData(rse_gene)

## read in lims
pheno = toxicant[[1]]
pheno = pheno[match(pd$BrNum, paste0("Br", pheno$brnumerical)),]
pd$BrainBank = ifelse(as.character(pheno$source)=="MD ME", "LIBD", "NIMH-HBCC")

# manifest
man = read.delim("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/.samples_unmerged.manifest",
	header=FALSE,as.is=TRUE)
man$RNum = ss(man$V5, "_")
man$BrNum = pd$BrNum[match(man$RNum, pd$RNum)]
man$BrainRegion = pd$BrainRegion[match(man$RNum, pd$RNum)]
man = man[,c(5,6,7,8,1,3)]

write.table(man, file = "fastqFiles_LIBD_bipolar.txt",
	sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	
## assay metadata file for RNA-seq
man$fileName = ss(sapply(strsplit(man$V1, "/"), function(x) x[length(x)]), "_R1")
pd$fileList = sapply(split(man$fileName, man$RNum)[pd$RNum], paste, collapse=";")

meta_rnaseq = data.frame(assay = "rnaSeq", individualID = pd$BrNum,
	individualIdSource = pd$BrainBank, fileName = pd$fileList,
	specimenID = pd$RNum, organ = "brain", tissue = pd$BrainRegion,
	brodmannArea = ifelse(pd$BrainRegion == "sACC", "BA25","Amygdala"),
	cellType = "Homogenate", RIN = pd$RIN, platform = "HiSeq3000",
	libraryPrep = "rRNAdepletion", runType ="pairedEnd", isStranded="true",
	readLength=100)
write.csv(meta_rnaseq, row.names=FALSE,
	file="ZandiHyde_BipSeq_LimbicData_RNAmetadata.csv")
	
## clinical data for BrNum
meta_brain = data.frame(individualID=pd$BrNum, age = pd$AgeDeath,
	sex = pd$Sex, race = pd$Race, ph = pd$pH, pmi = pd$PMI, 
	dx = pd$PrimaryDx, stringsAsFactors=FALSE)
meta_brain = meta_brain[!duplicated(meta_brain$individualID),]

write.csv(meta_brain, row.names=FALSE,
	file="ZandiHyde_BipSeq_LimbicData_SubjectLevel_metadata.csv")
	