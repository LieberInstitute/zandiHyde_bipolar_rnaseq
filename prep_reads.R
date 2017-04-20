###
library(jaffelab)

pd = read.csv("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/Bipolar_Zandi/Zandi RNA 20160323.csv")

## read path 1
theReads1 = list.files("/dcl01/lieber/ajaffe/Nina/Zandi/data", 
	pattern="fastq.gz$", full.names=TRUE)
theReads2 = list.files("/dcl01/ajaffe/data/Zandi_new/data", 
	pattern="fastq.gz$", full.names=TRUE)
theReads = c(theReads1, theReads2)	

## combine	
dat = data.frame(leftRead = theReads[grep("_R1_", theReads)], leftMd5=0,
	rightRead = theReads[grep("_R2_", theReads)], rightMd5=0,
	stringsAsFactors=FALSE)
	
## labels
dat$google = grepl("google", dat$leftRead)
dat$lab = gsub("/dcl01/ajaffe/data/Zandi_new/data/","", dat$leftRead)
dat$lab = gsub("/dcl01/lieber/ajaffe/Nina/Zandi/data/","", dat$lab)
dat$RNum = ss(dat$lab, "_")
dat$Flowcell = ss(dat$lab, "_",2)
dat$lab = paste0(dat$RNum, "_", dat$Flowcell)

dat = dat[,1:5]
write.table(dat, "preprocessed_data/samples.manifest", sep="\t",
	row.names=FALSE, col.names=FALSE, quote=FALSE)