load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/annotated_phenotype_data_zandi_hyde_bipolar_n511.rda")
all <- read.table("/users/schadinh/lieber/jxn_counts/BP_RNAseq_perind.counts.chrom.gz.qqnorm_all",header=T,comment.char="?")
cols <- colnames(all)
cols <- gsub("_junctions_primaryOnly_regtools","",cols)
br_all <- as.character(lapply(cols, function(x) ifelse(x %in% pd$SAMPLE_ID, pd[pd$SAMPLE_ID == x,"SampleID"],x)))
colnames(all) <- br_all
br_amy <- c(1, 2, 3, 4, grep("Amygdala",br_all))
amy <- all[,br_amy]
colnames(amy) <- gsub("_Amygdala","",colnames(amy))
write.table(amy, file="/users/schadinh/lieber/qqnorm/amygdala_qqnorm_20190329.txt", quote=F, sep="\t", col.names=T, row.names=F)
br_sacc <- c(1, 2, 3, 4, grep("sACC",br_all))
sacc <- all[,br_sacc]
colnames(sacc) <- gsub("_sACC","",colnames(sacc))
write.table(sacc, file="/users/schadinh/lieber/qqnorm/sacc_qqnorm_20190329.txt", quote=F, sep="\t", col.names=T, row.names=F)
