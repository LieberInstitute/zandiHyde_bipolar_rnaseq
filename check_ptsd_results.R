##

out = read.csv("ptsd_geneStats.csv", row.names=1,as.is=TRUE)
out$exprs_0.1_any = out$Amy_meanRpkm > 0.1 & out$Sacc_meanRpkm > 0.1 & out$DG_meanRpkm > 0.1

## meta analysis
out$meta_tstat_PTSD_CNT = rowSums(out[,grep("tstat_PTSD_CNT", names(out))])/sqrt(3)
out$meta_pvalue_PTSD_CNT = 2*pnorm(-abs(out$meta_tstat_PTSD_CNT ))
out$meta_bonf_PTSD_CNT = p.adjust(out$meta_pvalue_PTSD_CNT, "bonf")

out$meta_tstat_PTSD_BP = rowSums(out[,grep("tstat_PTSD_BP", names(out))])/sqrt(3)
out$meta_pvalue_PTSD_BP = 2*pnorm(-abs(out$meta_tstat_PTSD_BP ))
out$meta_bonf_PTSD_BP = p.adjust(out$meta_pvalue_PTSD_BP, "bonf")

## plot
plot(out$meta_tstat_PTSD_CNT ~ out$meta_tstat_PTSD_BP, 
	subset = out$exprs_0.1_any )

out$numMarg_PTSD_CNT = rowSums(out[,grep("pvalue_PTSD_CNT",names(out))[1:3]] < 1e-4)
out$numMarg_PTSD_BPD = rowSums(out[,grep("pvalue_PTSD_BP",names(out))[1:3]] < 1e-4)

## by region
sig = out[which(out$meta_bonf_PTSD_CNT < 0.2 | out$numMarg_PTSD_CNT > 0),]

## filter columns
sig = sig[,c(41, 44, 46:47, 6:10, 17:20, 25:28, 33:36, 12:14)]
sig = sig[order(sig$meta_pvalue_PTSD_CNT),]
write.csv(sig, file = "ptsd_geneStats_filtered.csv")


sig["ENSG00000100505.13",]

sacc = out[out$Sacc_pvalue_PTSD_CNT < 1e-4,]

numSig = rowSums(abs(out[,grep("tstat_PTSD_CNT", names(out))]) > 3) 
sig = out[numSig==2,]

plot(sig[,grep("tstat_PTSD_CNT", names(sig))])