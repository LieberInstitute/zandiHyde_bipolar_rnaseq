#Something I've noticed: the pSI function that calculates the genes falling into each threshold (specificity index) changes slightly each time you run it...



#Thanks for your help with this...can you run CSEA with the same 3 reference groups as those PTSD analyses (bacTRAP, then two human sets) for:
 
##DE genes: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda – you can stratify by features and do 
	#a similar sensitivity analysis as before to figure out what significance cutoff to use to get enough genes for the CSEA


#DE Genes
######Ask if want comparable cutoff values for the same feature...ie genes from sacc and amyg both use same cutoff...I prioritized getting same number genes but can prioritize same cutoff

library(jaffelab)
library(SummarizedExperiment)
library(pSI)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda")

#Only keep class = InGen and symbol != blank so we don't have issues later with selecting exactly 225 genes
statOut <- statOut[statOut$Class == "InGen" & statOut$Symbol != "",]

#Subset the features
statOutGene <- statOut[statOut$Type == "Gene",]
statOutExon <- statOut[statOut$Type == "Exon",]
statOutJxn <- statOut[statOut$Type == "Junction",]
statOutTx <- statOut[statOut$Type == "Transcript",]
dim(statOutGene)
#[1] 19205    17
dim(statOutExon)
#[1] 382411     17
dim(statOutJxn)
#[1] 208163     17
dim(statOutTx) 
#[1] 73214    17

#################################### Not Used #####################################

#Sensitivity analysis to find optimal cutoff for CSEA, Amyg
dim(statOutGene[statOutGene$adj.P.Val_Amyg < 0.07, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutExon[statOutExon$adj.P.Val_Amyg < 0.05, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutJxn[statOutJxn$adj.P.Val_Amyg < 0.25, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutTx[statOutTx$adj.P.Val_Amyg < 0.33, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")])


#Sensitivity analysis to find optimal cutoff for CSEA, sACC
dim(statOutGene[statOutGene$adj.P.Val_sACC < 0.014, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutExon[statOutExon$adj.P.Val_sACC < 0.008, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutJxn[statOutJxn$adj.P.Val_sACC < 0.12, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")])
dim(statOutTx[statOutTx$adj.P.Val_sACC < 0.09, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")])

#Name objects
gene_amyg_q07 <- statOutGene[statOutGene$adj.P.Val_Amyg < 0.07, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
exon_amyg_q05 <- statOutExon[statOutExon$adj.P.Val_Amyg < 0.05, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_amyg_q25 <- statOutJxn[statOutJxn$adj.P.Val_Amyg < 0.25, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
tx_amyg_q33 <- statOutTx[statOutTx$adj.P.Val_Amyg < 0.33, c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]

gene_sacc_q014 <- statOutGene[statOutGene$adj.P.Val_sACC < 0.014, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
exon_sacc_q008 <- statOutExon[statOutExon$adj.P.Val_sACC < 0.008, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_sacc_q12 <- statOutJxn[statOutJxn$adj.P.Val_sACC < 0.12, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
tx_sacc_q09 <- statOutTx[statOutTx$adj.P.Val_sACC < 0.09, c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]

#Save objects as rdas and csvs
save(gene_amyg_q07, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_q07.rda")
save(exon_amyg_q05, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_q05.rda")
save(jxn_amyg_q25, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_q25.rda")
save(tx_amyg_q33, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_q33.rda")

save(gene_sacc_q014, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_q014.rda")
save(exon_sacc_q008, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_q008.rda")
save(jxn_sacc_q12, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_q12.rda")
save(tx_sacc_q09, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_q09.rda")

write.csv(gene_amyg_q07, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_q07.csv")
write.csv(exon_amyg_q05, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_Amyg_q05.csv")
write.csv(jxn_amyg_q25, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_Amyg_q25.csv")
write.csv(tx_amyg_q33, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_Amyg_q33.csv")

write.csv(gene_sacc_q014, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_q014.csv")
write.csv(exon_sacc_q008, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_sACC_q008.csv")
write.csv(jxn_sacc_q12, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_sACC_q12.csv")
write.csv(tx_sacc_q09, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_sACC_q09.csv")


##############################################################################################################################

#Well, the above doesn't count only the unique genes. Proper sensitivity analysis would look like this:
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.08, c("logFC_sACC","Symbol")]$Symbol))


###############################################################################################
########################Sensitivity analysis to see best cutoff################################
###############################################################################################
#sACC genes

length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.5, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.4, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.3, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.2, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.1, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.05, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.01, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.005, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.001, "Symbol"]))

genes_sACC_q5 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.5, "Symbol"])
genes_sACC_q4 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.4, "Symbol"])
genes_sACC_q3 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.3, "Symbol"])
genes_sACC_q2 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.2, "Symbol"])
genes_sACC_q1 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.1, "Symbol"])
genes_sACC_q05 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.05, "Symbol"])
genes_sACC_q01 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.01, "Symbol"])
genes_sACC_q005 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.005, "Symbol"])
genes_sACC_q001 <- unique(statOutGene[statOutGene$adj.P.Val_sACC < 0.001, "Symbol"])

#sACC exons

length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.4, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.3, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.2, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.1, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.05, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.01, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.005, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.001, "Symbol"]))

exons_sACC_q4 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.4, "Symbol"])
exons_sACC_q3 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.3, "Symbol"])
exons_sACC_q2 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.2, "Symbol"])
exons_sACC_q1 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.1, "Symbol"])
exons_sACC_q05 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.05, "Symbol"])
exons_sACC_q01 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.01, "Symbol"])
exons_sACC_q005 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.005, "Symbol"])
exons_sACC_q001 <- unique(statOutExon[statOutExon$adj.P.Val_sACC < 0.001, "Symbol"])




#Amyg genes (not as many tests bc not as many genes as sACC)
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.5, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.4, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.3, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.2, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.1, "Symbol"]))
length(unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.05, "Symbol"]))

genes_Amyg_q5 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.5, "Symbol"])
genes_Amyg_q4 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.4, "Symbol"])
genes_Amyg_q3 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.3, "Symbol"])
genes_Amyg_q2 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.2, "Symbol"])
genes_Amyg_q1 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.1, "Symbol"])
genes_Amyg_q05 <- unique(statOutGene[statOutGene$adj.P.Val_Amyg < 0.05, "Symbol"])

#Amyg exons
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.5, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.4, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.3, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.2, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.1, "Symbol"]))
length(unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.05, "Symbol"]))

exons_Amyg_q5 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.5, "Symbol"])
exons_Amyg_q4 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.4, "Symbol"])
exons_Amyg_q3 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.3, "Symbol"])
exons_Amyg_q2 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.2, "Symbol"])
exons_Amyg_q1 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.1, "Symbol"])
exons_Amyg_q05 <- unique(statOutExon[statOutExon$adj.P.Val_Amyg < 0.05, "Symbol"])

##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined
##Broader cell type definitions, still doing across regions

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg_broadtypes.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed_broad[is.na(avg_lake_2018_all_renamed_broad)] <- 0
avg_lake_2018_all_renamed_broad_notfiltered <- avg_lake_2018_all_renamed_broad
avg_lake_2018_all_renamed_broad[avg_lake_2018_all_renamed_broad < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_broad_notfiltered, avg_lake_2018_all_renamed_broad, p_max=1)

pSI.count(lake_spi)	

##Calculation

#sACC genes
fisher.iteration(lake_spi,genes_sACC_q5)
fisher.iteration(lake_spi,genes_sACC_q4)
fisher.iteration(lake_spi,genes_sACC_q3)
fisher.iteration(lake_spi,genes_sACC_q2)
fisher.iteration(lake_spi,genes_sACC_q1)
fisher.iteration(lake_spi,genes_sACC_q05)
fisher.iteration(lake_spi,genes_sACC_q01)
fisher.iteration(lake_spi,genes_sACC_q005)
fisher.iteration(lake_spi,genes_sACC_q001)

#sACC exons
fisher.iteration(lake_spi,exons_sACC_q4)
fisher.iteration(lake_spi,exons_sACC_q3)
fisher.iteration(lake_spi,exons_sACC_q2)
fisher.iteration(lake_spi,exons_sACC_q1)
fisher.iteration(lake_spi,exons_sACC_q05)
fisher.iteration(lake_spi,exons_sACC_q01)
fisher.iteration(lake_spi,exons_sACC_q005)
fisher.iteration(lake_spi,exons_sACC_q001)


#Amyg genes
fisher.iteration(lake_spi,genes_Amyg_q5)
fisher.iteration(lake_spi,genes_Amyg_q4)
fisher.iteration(lake_spi,genes_Amyg_q3)
fisher.iteration(lake_spi,genes_Amyg_q2)
fisher.iteration(lake_spi,genes_Amyg_q1)
fisher.iteration(lake_spi,genes_Amyg_q05)

#Amyg exons
fisher.iteration(lake_spi,exons_Amyg_q5)
fisher.iteration(lake_spi,exons_Amyg_q4)
fisher.iteration(lake_spi,exons_Amyg_q3)
fisher.iteration(lake_spi,exons_Amyg_q2)
fisher.iteration(lake_spi,exons_Amyg_q1)
fisher.iteration(lake_spi,exons_Amyg_q05)

##Save output/results

#sACC genes
write.csv(fisher.iteration(lake_spi,genes_sACC_q5), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q5_Lake_broad.csv")
write.csv(fisher.iteration(lake_spi,genes_sACC_q4), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q4_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q3), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q3_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q2), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q2_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q1), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q1_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q05), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q05_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q01), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q01_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q005), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q005_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_sACC_q001), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_sACC_q001_Lake_broad.csv")				

#sACC exons
write.csv(fisher.iteration(lake_spi,exons_sACC_q4), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q4_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q3), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q3_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q2), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q2_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q1), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q1_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q05), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q05_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q01), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q01_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q005), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q005_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_sACC_q001), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_sACC_q001_Lake_broad.csv")				


#Amyg genes
write.csv(fisher.iteration(lake_spi,genes_Amyg_q5), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q5_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_Amyg_q4), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q4_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_Amyg_q3), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q3_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_Amyg_q2), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q2_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_Amyg_q1), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q1_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_Amyg_q05), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/genes_Amyg_q05_Lake_broad.csv")				

#Amyg exons
write.csv(fisher.iteration(lake_spi,exons_Amyg_q5), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q5_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_Amyg_q4), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q4_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_Amyg_q3), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q3_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_Amyg_q2), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q2_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_Amyg_q1), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q1_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exons_Amyg_q05), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/Sensitivity_Analyses/exons_Amyg_q05_Lake_broad.csv")				


##############################################################################################################################


#I'm effectively just trying to take the top ~225 genes, and since the q and p values differ so much between regions in the number of genes, let's just take the top 225 unique genes based on ordered q values
#Note, some of these features will have a gencodeID, etc. but no gene symbol, and since CSEA takes gene symbols only, I took the top 225 unique gene symbol names


#Order objects based on q values
statOutGeneAmyg <- statOutGene[order(statOutGene$adj.P.Val_Amyg),]
statOutExonAmyg <- statOutExon[order(statOutExon$adj.P.Val_Amyg),]
statOutJxnAmyg <- statOutJxn[order(statOutJxn$adj.P.Val_Amyg),]
statOutTxAmyg <- statOutTx[order(statOutTx$adj.P.Val_Amyg),]

statOutGenesACC <- statOutGene[order(statOutGene$adj.P.Val_sACC),]
statOutExonsACC <- statOutExon[order(statOutExon$adj.P.Val_sACC),]
statOutJxnsACC <- statOutJxn[order(statOutJxn$adj.P.Val_sACC),]
statOutTxsACC <- statOutTx[order(statOutTx$adj.P.Val_sACC),]

#find number of rows needed for 225 unique genes 
length(unique(statOutGeneAmyg[c(1:225),"Symbol"]))
length(unique(statOutExonAmyg[c(1:577),"Symbol"]))
length(unique(statOutJxnAmyg[c(1:312),"Symbol"]))
length(unique(statOutTxAmyg[c(1:239),"Symbol"]))

length(unique(statOutGenesACC[c(1:225),"Symbol"]))
length(unique(statOutExonsACC[c(1:611),"Symbol"]))
length(unique(statOutJxnsACC[c(1:362),"Symbol"]))
length(unique(statOutTxsACC[c(1:237),"Symbol"]))


#Name objects
gene_amyg_225 <- statOutGeneAmyg[c(1:225), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
exon_amyg_225 <- statOutExonAmyg[c(1:577), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_amyg_225 <- statOutJxnAmyg[c(1:312), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]
tx_amyg_225 <- statOutTxAmyg[c(1:239), c("logFC_Amyg","AveExpr_Amyg","t_Amyg","P.Value_Amyg","adj.P.Val_Amyg","B_Amyg","gencodeID","Symbol","Type","Class","EntrezID")]

gene_sacc_225 <- statOutGenesACC[c(1:225), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
exon_sacc_225 <- statOutExonsACC[c(1:611), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
jxn_sacc_225 <- statOutJxnsACC[c(1:362), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]
tx_sacc_225 <- statOutTxsACC[c(1:237), c("logFC_sACC","AveExpr_sACC","t_sACC","P.Value_sACC","adj.P.Val_sACC","B_sACC","gencodeID","Symbol","Type","Class","EntrezID")]

#Save objects as rdas and csvs
save(gene_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
save(exon_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
save(jxn_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
save(tx_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")

save(gene_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
save(exon_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
save(jxn_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
save(tx_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")

write.csv(gene_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_225.csv")
write.csv(exon_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_Amyg_225.csv")
write.csv(jxn_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_Amyg_225.csv")
write.csv(tx_amyg_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_Amyg_225.csv")

write.csv(gene_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_225.csv")
write.csv(exon_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_sACC_225.csv")
write.csv(jxn_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_sACC_225.csv")
write.csv(tx_sacc_225, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_sACC_225.csv")




##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### CSEA
### Online tool (http://genetics.wustl.edu/jdlab/csea-tool-2/)
##Use csvs generated above

#Note, CSEA may not be the best for things other than genes because each gene symbol is only counted once...what if multiple transcripts map to the same gene? Shouldn't it count more?



##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.
##Need to select more specific p value?
##need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done... filter before or after average calculation?
##need to do any normalization first?
##UMI?

library(jaffelab)
library(SummarizedExperiment)
library(pSI)


#get data for pSI calculation
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Li_2018_renamed.RData")

#need to filter out lowly expressed stuff?
named_not_filtered <- new2_umi2_renamed
new2_umi2_renamed[new2_umi2_renamed < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
hm <- specificity.index(named_not_filtered, new2_umi2_renamed, p_max=1)

pSI.count(hm)
pSI.list(hm, write.csv=FALSE)		


#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")

gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(hm,gene_amyg_225_sym)				
fisher.iteration(hm,exon_amyg_225_sym)
fisher.iteration(hm,jxn_amyg_225_sym)	
fisher.iteration(hm,tx_amyg_225_sym)		

write.csv(fisher.iteration(hm,gene_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Li.csv")				
write.csv(fisher.iteration(hm,exon_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Li.csv")		
write.csv(fisher.iteration(hm,jxn_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Li.csv")		
write.csv(fisher.iteration(hm,tx_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Li.csv")		


#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")

gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(hm,gene_sacc_225_sym)			
fisher.iteration(hm,exon_sacc_225_sym)		
fisher.iteration(hm,jxn_sacc_225_sym)		
fisher.iteration(hm,tx_sacc_225_sym)		

write.csv(fisher.iteration(hm,gene_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Li.csv")				
write.csv(fisher.iteration(hm,exon_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Li.csv")		
write.csv(fisher.iteration(hm,jxn_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Li.csv")		
write.csv(fisher.iteration(hm,tx_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Li.csv")		
				




##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed[is.na(avg_lake_2018_all_renamed)] <- 0
avg_lake_2018_all_renamed_notfiltered <- avg_lake_2018_all_renamed
avg_lake_2018_all_renamed[avg_lake_2018_all_renamed < .2] <- NA


#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_notfiltered, avg_lake_2018_all_renamed, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		

#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")

gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(lake_spi,gene_amyg_225_sym)			
fisher.iteration(lake_spi,exon_amyg_225_sym)	
fisher.iteration(lake_spi,jxn_amyg_225_sym)		
fisher.iteration(lake_spi,tx_amyg_225_sym)		

write.csv(fisher.iteration(lake_spi,gene_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Lake.csv")				
write.csv(fisher.iteration(lake_spi,exon_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Lake.csv")		
write.csv(fisher.iteration(lake_spi,jxn_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Lake.csv")		
write.csv(fisher.iteration(lake_spi,tx_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Lake.csv")		

#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")

gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(lake_spi,gene_sacc_225_sym)				
fisher.iteration(lake_spi,exon_sacc_225_sym)		
fisher.iteration(lake_spi,jxn_sacc_225_sym)		
fisher.iteration(lake_spi,tx_sacc_225_sym)		


write.csv(fisher.iteration(lake_spi,gene_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Lake.csv")				
write.csv(fisher.iteration(lake_spi,exon_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Lake.csv")		
write.csv(fisher.iteration(lake_spi,jxn_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Lake.csv")		
write.csv(fisher.iteration(lake_spi,tx_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Lake.csv")		

##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined
##Broader cell type definitions, still doing across regions

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg_broadtypes.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed_broad[is.na(avg_lake_2018_all_renamed_broad)] <- 0
avg_lake_2018_all_renamed_broad_notfiltered <- avg_lake_2018_all_renamed_broad
avg_lake_2018_all_renamed_broad[avg_lake_2018_all_renamed_broad < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_broad_notfiltered, avg_lake_2018_all_renamed_broad, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		


#Amyg, all features

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_225.rda")

gene_amyg_225_sym <- gene_amyg_225[,"Symbol"]
exon_amyg_225_sym <- exon_amyg_225[,"Symbol"]
jxn_amyg_225_sym <- jxn_amyg_225[,"Symbol"]
tx_amyg_225_sym <- tx_amyg_225[,"Symbol"]

fisher.iteration(lake_spi,gene_amyg_225_sym)			
fisher.iteration(lake_spi,exon_amyg_225_sym)	
fisher.iteration(lake_spi,jxn_amyg_225_sym)		
fisher.iteration(lake_spi,tx_amyg_225_sym)		

write.csv(fisher.iteration(lake_spi,gene_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_225_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exon_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_225_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,jxn_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_225_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,tx_amyg_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_225_Lake_broad.csv")		

#sACC

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_225.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_225.rda")

gene_sacc_225_sym <- gene_sacc_225[,"Symbol"]
exon_sacc_225_sym <- exon_sacc_225[,"Symbol"]
jxn_sacc_225_sym <- jxn_sacc_225[,"Symbol"]
tx_sacc_225_sym <- tx_sacc_225[,"Symbol"]

fisher.iteration(lake_spi,gene_sacc_225_sym)				
fisher.iteration(lake_spi,exon_sacc_225_sym)		
fisher.iteration(lake_spi,jxn_sacc_225_sym)		
fisher.iteration(lake_spi,tx_sacc_225_sym)		


write.csv(fisher.iteration(lake_spi,gene_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_225_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exon_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_225_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,jxn_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_225_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,tx_sacc_225_sym), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_225_Lake_broad.csv")		



##eQTL genes: you can follow leo’s script for GO here: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr/go_analysis.R to see what genes went in, but lets just 
	#start with the genes significantly enumerated on lines 105-106 (263/250/381)


#eQTL genes
#using gene symbols since that's the only input for CSEA

library('SummarizedExperiment')
library('data.table')
library('clusterProfiler')
library('GenomicRanges')
library('sessioninfo')
library('jaffelab')
library('pSI')

setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl/raggr")


## Load the 881 eQTL results and get the unique gene IDs
eqtl_files <- dir(pattern = 'raggr_881_snps_.*')
names(eqtl_files) <- gsub('raggr_881_snps_|_eqtls_fdr01.csv', '', eqtl_files)

#subset into the different features
genes <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Gene",]
})

exons <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Exon",]
})

junctions <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Jxn",]
})

transcripts <- lapply(eqtl_files, function(ef) {
	qtl <- fread(ef)
	qtl[qtl$Type == "Tx",]
})


#get the unique, non-blank gene symbols for each feature

sig_genes <- lapply(genes,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_genes, length)
# amyg dlpfc  sacc 
#   84    89   144 

sig_exons <- lapply(exons,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_exons, length)
#amyg dlpfc  sacc 
#  132   126   208 

sig_junctions <- lapply(junctions,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_junctions, length)
#amyg dlpfc  sacc 
#  119   106   160 
sig_transcripts <- lapply(transcripts,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_transcripts, length)
#amyg dlpfc  sacc 
#   90   121   142 


#Subset the gene symbols by region

gene_amyg_eqtl <- sig_genes$amyg
exon_amyg_eqtl <- sig_exons$amyg
jxn_amyg_eqtl <- sig_junctions$amyg
tx_amyg_eqtl <- sig_transcripts$amyg

gene_dlpfc_eqtl <- sig_genes$dlpfc
exon_dlpfc_eqtl <- sig_exons$dlpfc
jxn_dlpfc_eqtl <- sig_junctions$dlpfc
tx_dlpfc_eqtl <- sig_transcripts$dlpfc

gene_sacc_eqtl <- sig_genes$sacc
exon_sacc_eqtl <- sig_exons$sacc
jxn_sacc_eqtl <- sig_junctions$sacc
tx_sacc_eqtl <- sig_transcripts$sacc


#Save objects as rdas and csvs
save(gene_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
save(exon_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
save(jxn_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
save(tx_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")

save(gene_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
save(exon_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
save(jxn_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
save(tx_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")

save(gene_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
save(exon_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
save(jxn_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
save(tx_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")

write.csv(gene_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_eQTL.csv")
write.csv(exon_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_Amyg_eQTL.csv")
write.csv(jxn_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_Amyg_eQTL.csv")
write.csv(tx_amyg_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_Amyg_eQTL.csv")

write.csv(gene_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_DLPFC_eQTL.csv")
write.csv(exon_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_DLPFC_eQTL.csv")
write.csv(jxn_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_DLPFC_eQTL.csv")
write.csv(tx_dlpfc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_DLPFC_eQTL.csv")

write.csv(gene_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_eQTL.csv")
write.csv(exon_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/exons_Symbols_sACC_eQTL.csv")
write.csv(jxn_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/junctions_Symbols_sACC_eQTL.csv")
write.csv(tx_sacc_eqtl, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/transcripts_Symbols_sACC_eQTL.csv")



#all genes from all features together

## Load the 881 eQTL results and get the unique, not blank gene symbols from all features
eqtl_files <- dir(pattern = 'raggr_881_snps_.*')
names(eqtl_files) <- gsub('raggr_881_snps_|_eqtls_fdr01.csv', '', eqtl_files)

qtl <- lapply(eqtl_files, fread)
})

sig_genes_all <- lapply(qtl,function(x) {
	x$Symbol[x$Symbol == ""] <- NA
	unique(x$Symbol[!is.na(x$Symbol)])
})
sapply(sig_genes_all, length)
# amyg dlpfc  sacc 
#  228   249   328 

gene_amyg_eqtl_all <- sig_genes_all$amyg
gene_dlpfc_eqtl_all <- sig_genes_all$dlpfc
gene_sacc_eqtl_all <- sig_genes_all$sacc


save(gene_amyg_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")
save(gene_dlpfc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")
save(gene_sacc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

write.csv(gene_amyg_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_Amyg_eQTL_allFeatures.csv")
write.csv(gene_dlpfc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_DLPFC_eQTL_allFeatures.csv")
write.csv(gene_sacc_eqtl_all, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/genes_Symbols_sACC_eQTL_allFeatures.csv")


##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### CSEA
### Online tool (http://genetics.wustl.edu/jdlab/csea-tool-2/)
##Use csvs generated above

#Note, CSEA may not be the best for things other than genes because each gene symbol is only counted once...what if multiple transcripts map to the same gene? Shouldn't it count more?




##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.
##Need to select more specific p value?
##need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done... filter before or after average calculation?
##need to do any normalization first?
##UMI?

library(jaffelab)
library(SummarizedExperiment)
library(pSI)


#get data for pSI calculation
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Li_2018_renamed.RData")

#need to filter out lowly expressed stuff?
named_not_filtered <- new2_umi2_renamed
new2_umi2_renamed[new2_umi2_renamed < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
hm <- specificity.index(named_not_filtered, new2_umi2_renamed, p_max=1)

pSI.count(hm)
pSI.list(hm, write.csv=FALSE)		


#Amyg, all features 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")

fisher.iteration(hm,gene_amyg_eqtl)				
fisher.iteration(hm,exon_amyg_eqtl)
fisher.iteration(hm,jxn_amyg_eqtl)	
fisher.iteration(hm,tx_amyg_eqtl)		
fisher.iteration(hm,gene_amyg_eqtl_all)		

write.csv(fisher.iteration(hm,gene_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Li.csv")				
write.csv(fisher.iteration(hm,exon_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,jxn_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,tx_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,gene_amyg_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Li.csv")		


#DLPFC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")

fisher.iteration(hm,gene_dlpfc_eqtl)				
fisher.iteration(hm,exon_dlpfc_eqtl)
fisher.iteration(hm,jxn_dlpfc_eqtl)	
fisher.iteration(hm,tx_dlpfc_eqtl)		
fisher.iteration(hm,gene_dlpfc_eqtl_all)		

write.csv(fisher.iteration(hm,gene_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Li.csv")				
write.csv(fisher.iteration(hm,exon_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,jxn_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,tx_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,gene_dlpfc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Li.csv")		


#sACC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

fisher.iteration(hm,gene_sacc_eqtl)				
fisher.iteration(hm,exon_sacc_eqtl)
fisher.iteration(hm,jxn_sacc_eqtl)	
fisher.iteration(hm,tx_sacc_eqtl)		
fisher.iteration(hm,gene_sacc_eqtl_all)		

write.csv(fisher.iteration(hm,gene_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Li.csv")				
write.csv(fisher.iteration(hm,exon_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,jxn_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,tx_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Li.csv")		
write.csv(fisher.iteration(hm,gene_sacc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Li.csv")		


##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed[is.na(avg_lake_2018_all_renamed)] <- 0
avg_lake_2018_all_renamed_notfiltered <- avg_lake_2018_all_renamed
avg_lake_2018_all_renamed[avg_lake_2018_all_renamed < .2] <- NA


#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_notfiltered, avg_lake_2018_all_renamed, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		

#Amyg, all features 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_amyg_eqtl)				
fisher.iteration(lake_spi,exon_amyg_eqtl)
fisher.iteration(lake_spi,jxn_amyg_eqtl)	
fisher.iteration(lake_spi,tx_amyg_eqtl)		
fisher.iteration(lake_spi,gene_amyg_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Lake.csv")				
write.csv(fisher.iteration(lake_spi,exon_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,jxn_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,tx_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,gene_amyg_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Lake.csv")		


#DLPFC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_dlpfc_eqtl)				
fisher.iteration(lake_spi,exon_dlpfc_eqtl)
fisher.iteration(lake_spi,jxn_dlpfc_eqtl)	
fisher.iteration(lake_spi,tx_dlpfc_eqtl)		
fisher.iteration(lake_spi,gene_dlpfc_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Lake.csv")				
write.csv(fisher.iteration(lake_spi,exon_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,jxn_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,tx_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,gene_dlpfc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Lake.csv")		


#sACC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_sacc_eqtl)				
fisher.iteration(lake_spi,exon_sacc_eqtl)
fisher.iteration(lake_spi,jxn_sacc_eqtl)	
fisher.iteration(lake_spi,tx_sacc_eqtl)		
fisher.iteration(lake_spi,gene_sacc_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Lake.csv")				
write.csv(fisher.iteration(lake_spi,exon_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,jxn_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,tx_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Lake.csv")		
write.csv(fisher.iteration(lake_spi,gene_sacc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Lake.csv")		



##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined
##Broader cell type definitions, still doing across regions

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg_broadtypes.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed_broad[is.na(avg_lake_2018_all_renamed_broad)] <- 0
avg_lake_2018_all_renamed_broad_notfiltered <- avg_lake_2018_all_renamed_broad
avg_lake_2018_all_renamed_broad[avg_lake_2018_all_renamed_broad < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_broad_notfiltered, avg_lake_2018_all_renamed_broad, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		

#Amyg, all features 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_Amyg_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_Amyg_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_amyg_eqtl)				
fisher.iteration(lake_spi,exon_amyg_eqtl)
fisher.iteration(lake_spi,jxn_amyg_eqtl)	
fisher.iteration(lake_spi,tx_amyg_eqtl)		
fisher.iteration(lake_spi,gene_amyg_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exon_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_Amyg_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,jxn_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_Amyg_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,tx_amyg_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_Amyg_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,gene_amyg_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_Amyg_eQTL_allFeatures_Lake_broad.csv")		


#DLPFC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_DLPFC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_DLPFC_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_dlpfc_eqtl)				
fisher.iteration(lake_spi,exon_dlpfc_eqtl)
fisher.iteration(lake_spi,jxn_dlpfc_eqtl)	
fisher.iteration(lake_spi,tx_dlpfc_eqtl)		
fisher.iteration(lake_spi,gene_dlpfc_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exon_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_DLPFC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,jxn_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_DLPFC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,tx_dlpfc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_DLPFC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,gene_dlpfc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_DLPFC_eQTL_allFeatures_Lake_broad.csv")		


#sACC, all features
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/exons_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/junctions_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/transcripts_Symbols_sACC_eQTL.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/genes_Symbols_sACC_eQTL_allFeatures.rda")

fisher.iteration(lake_spi,gene_sacc_eqtl)				
fisher.iteration(lake_spi,exon_sacc_eqtl)
fisher.iteration(lake_spi,jxn_sacc_eqtl)	
fisher.iteration(lake_spi,tx_sacc_eqtl)		
fisher.iteration(lake_spi,gene_sacc_eqtl_all)		

write.csv(fisher.iteration(lake_spi,gene_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,exon_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/exons_sACC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,jxn_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/junctions_sACC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,tx_sacc_eqtl), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/transcripts_sACC_eQTL_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,gene_sacc_eqtl_all), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_sACC_eQTL_allFeatures_Lake_broad.csv")		



##WGCNA modules: just which genes are in which modules to assign a cell type enrichment to each module: /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/wgcna/analyze_results.R 
	#you can see how I got the genes for GO enrichment on lines 65-67

#this analysis is across brain regions (amygdala and sACC)
#this might be too many genes...for example, grey module has unique, non-blank 8817 genes
#other modules may have too few (~55)

library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
library(clusterProfiler)
library(readxl)

setwd("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/")

## load data
load("data/zandiHypde_bipolar_rseGene_n511.rda")
load("data/degradation_rse_BipSeq_BothRegions.rda")

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"),
                                levels = c("Control", "Bipolar"))

## add ancestry
load("genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

###########
# filter ##
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]
#number of gene symbols: 25136
#unique = 19201
#unique and non-blank gene symbols = 19200


#########################
## load wgcna output ####
#########################

## load
load("wgcna/rdas/constructed_network_signed_bicor.rda")

# get colors
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab,
        stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

#convert blank gene symbols to NA
length(rowData(rse_gene)$Symbol)
#25136
rowData(rse_gene)$Symbol[rowData(rse_gene)$Symbol == ""] <- NA
length(rowData(rse_gene)$Symbol)
#25136

#get unique, non-NA gene symbols
gList = split(rowData(rse_gene)$Symbol,
        factor(net$colorsLab,levels = colorDat$col))
sapply(gList, length)
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
sapply(gList, length)
gList = lapply(gList, function(x) as.character(unique(x[!is.na(x)])))
sapply(gList, length)
#the unique function only changed the grey module slightly, rest stayed same

#		 grey    turquoise         blue        brown       yellow        green 
#        8817         2147         2191         1489         1215         1043 
#         red        black         pink      magenta       purple  greenyellow 
#         518          273          216          217          196          180 
#         tan       salmon         cyan midnightblue    lightcyan       grey60 
#         137           93           75           75           77           59 
#  lightgreen  lightyellow    royalblue 
#          61           65           57 



#separate each module

genes_wgcna_grey <- gList$grey
genes_wgcna_turquoise <- gList$turquoise
genes_wgcna_blue <- gList$blue
genes_wgcna_brown <- gList$brown
genes_wgcna_yellow <- gList$yellow
genes_wgcna_green <- gList$green
genes_wgcna_red <- gList$red
genes_wgcna_black <- gList$black
genes_wgcna_pink <- gList$pink
genes_wgcna_magenta <- gList$magenta
genes_wgcna_purple <- gList$purple
genes_wgcna_greenyellow <- gList$greenyellow
genes_wgcna_tan <- gList$tan
genes_wgcna_salmon <- gList$salmon
genes_wgcna_cyan <- gList$cyan
genes_wgcna_midnightblue <- gList$midnightblue
genes_wgcna_lightcyan <- gList$lightcyan
genes_wgcna_grey60 <- gList$grey60
genes_wgcna_lightgreen <- gList$lightgreen
genes_wgcna_lightyellow <- gList$lightyellow
genes_wgcna_royalblue <- gList$royalblue

#spotcheck
length(genes_wgcna_grey)
#[1] 8817
length(genes_wgcna_turquoise)
#[1] 2147
length(genes_wgcna_pink)
#[1] 216
length(genes_wgcna_magenta)
#[1] 217
length(genes_wgcna_lightcyan)
#[1] 77
length(genes_wgcna_royalblue)
#[1] 57


#Save objects as rdas and csvs

save(genes_wgcna_grey, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
save(genes_wgcna_turquoise, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
save(genes_wgcna_blue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
save(genes_wgcna_brown, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
save(genes_wgcna_yellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
save(genes_wgcna_green, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
save(genes_wgcna_red, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
save(genes_wgcna_black, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
save(genes_wgcna_pink, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
save(genes_wgcna_magenta, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
save(genes_wgcna_purple, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
save(genes_wgcna_greenyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
save(genes_wgcna_tan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
save(genes_wgcna_salmon, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
save(genes_wgcna_cyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
save(genes_wgcna_midnightblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
save(genes_wgcna_lightcyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
save(genes_wgcna_grey60, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
save(genes_wgcna_lightgreen, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
save(genes_wgcna_lightyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
save(genes_wgcna_royalblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

write.csv(genes_wgcna_grey, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_grey.csv")
write.csv(genes_wgcna_turquoise, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_turquoise.csv")
write.csv(genes_wgcna_blue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_blue.csv")
write.csv(genes_wgcna_brown, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_brown.csv")
write.csv(genes_wgcna_yellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_yellow.csv")
write.csv(genes_wgcna_green, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_green.csv")
write.csv(genes_wgcna_red, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_red.csv")
write.csv(genes_wgcna_black, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_black.csv")
write.csv(genes_wgcna_pink, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_pink.csv")
write.csv(genes_wgcna_magenta, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_magenta.csv")
write.csv(genes_wgcna_purple, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_purple.csv")
write.csv(genes_wgcna_greenyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_greenyellow.csv")
write.csv(genes_wgcna_tan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_tan.csv")
write.csv(genes_wgcna_salmon, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_salmon.csv")
write.csv(genes_wgcna_cyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_cyan.csv")
write.csv(genes_wgcna_midnightblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_midnightblue.csv")
write.csv(genes_wgcna_lightcyan, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightcyan.csv")
write.csv(genes_wgcna_grey60, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_grey60.csv")
write.csv(genes_wgcna_lightgreen, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightgreen.csv")
write.csv(genes_wgcna_lightyellow, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_lightyellow.csv")
write.csv(genes_wgcna_royalblue, file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/gene_Symbols_WGCNA_royalblue.csv")




##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### CSEA
### Online tool (http://genetics.wustl.edu/jdlab/csea-tool-2/)
##Use csvs generated above

#Note, CSEA may not be the best for things other than genes because each gene symbol is only counted once...what if multiple transcripts map to the same gene? Shouldn't it count more?




##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.
##Need to select more specific p value?
##need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done... filter before or after average calculation?
##need to do any normalization first?
##UMI?

library(jaffelab)
library(SummarizedExperiment)
library(pSI)


#get data for pSI calculation
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Li_2018_renamed.RData")

#need to filter out lowly expressed stuff?
named_not_filtered <- new2_umi2_renamed
new2_umi2_renamed[new2_umi2_renamed < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
hm <- specificity.index(named_not_filtered, new2_umi2_renamed, p_max=1)

pSI.count(hm)
pSI.list(hm, write.csv=FALSE)		

#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

#Calculation

fisher.iteration(hm,genes_wgcna_grey)
fisher.iteration(hm,genes_wgcna_turquoise)	
fisher.iteration(hm,genes_wgcna_blue)	
fisher.iteration(hm,genes_wgcna_brown)	
fisher.iteration(hm,genes_wgcna_yellow)	
fisher.iteration(hm,genes_wgcna_green)	
fisher.iteration(hm,genes_wgcna_red)	
fisher.iteration(hm,genes_wgcna_black)	
fisher.iteration(hm,genes_wgcna_pink)	
fisher.iteration(hm,genes_wgcna_magenta)	
fisher.iteration(hm,genes_wgcna_purple)	
fisher.iteration(hm,genes_wgcna_greenyellow)	
fisher.iteration(hm,genes_wgcna_tan)	
fisher.iteration(hm,genes_wgcna_salmon)	
fisher.iteration(hm,genes_wgcna_cyan)	
fisher.iteration(hm,genes_wgcna_midnightblue)	
fisher.iteration(hm,genes_wgcna_lightcyan)	
fisher.iteration(hm,genes_wgcna_grey60)	
fisher.iteration(hm,genes_wgcna_lightgreen)	
fisher.iteration(hm,genes_wgcna_lightyellow)	
fisher.iteration(hm,genes_wgcna_royalblue)	
	
#Save output/results

write.csv(fisher.iteration(hm,genes_wgcna_grey), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Li.csv")				
write.csv(fisher.iteration(hm,genes_wgcna_turquoise), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_blue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_brown), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_yellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_green), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_red), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_black), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_pink), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_magenta), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_purple), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_greenyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_tan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_salmon), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_cyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_midnightblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_lightcyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_grey60), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_lightgreen), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_lightyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Li.csv")		
write.csv(fisher.iteration(hm,genes_wgcna_royalblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Li.csv")		



##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed[is.na(avg_lake_2018_all_renamed)] <- 0
avg_lake_2018_all_renamed_notfiltered <- avg_lake_2018_all_renamed
avg_lake_2018_all_renamed[avg_lake_2018_all_renamed < .2] <- NA


#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_notfiltered, avg_lake_2018_all_renamed, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		


#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

#Calculation

fisher.iteration(lake_spi,genes_wgcna_grey)
fisher.iteration(lake_spi,genes_wgcna_turquoise)	
fisher.iteration(lake_spi,genes_wgcna_blue)	
fisher.iteration(lake_spi,genes_wgcna_brown)	
fisher.iteration(lake_spi,genes_wgcna_yellow)	
fisher.iteration(lake_spi,genes_wgcna_green)	
fisher.iteration(lake_spi,genes_wgcna_red)	
fisher.iteration(lake_spi,genes_wgcna_black)	
fisher.iteration(lake_spi,genes_wgcna_pink)	
fisher.iteration(lake_spi,genes_wgcna_magenta)	
fisher.iteration(lake_spi,genes_wgcna_purple)	
fisher.iteration(lake_spi,genes_wgcna_greenyellow)	
fisher.iteration(lake_spi,genes_wgcna_tan)	
fisher.iteration(lake_spi,genes_wgcna_salmon)	
fisher.iteration(lake_spi,genes_wgcna_cyan)	
fisher.iteration(lake_spi,genes_wgcna_midnightblue)	
fisher.iteration(lake_spi,genes_wgcna_lightcyan)	
fisher.iteration(lake_spi,genes_wgcna_grey60)	
fisher.iteration(lake_spi,genes_wgcna_lightgreen)	
fisher.iteration(lake_spi,genes_wgcna_lightyellow)	
fisher.iteration(lake_spi,genes_wgcna_royalblue)	
	
#Save output/results

write.csv(fisher.iteration(lake_spi,genes_wgcna_grey), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Lake.csv")				
write.csv(fisher.iteration(lake_spi,genes_wgcna_turquoise), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_blue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_brown), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_yellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_green), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_red), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_black), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_pink), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_magenta), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_purple), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_greenyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_tan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_salmon), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_cyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_midnightblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightcyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_grey60), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightgreen), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Lake.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_royalblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Lake.csv")		
	



##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers for frontal cotex, visual cortex, and cerebellum separately and combined
##Broader cell type definitions, still doing across regions

library(jaffelab)
library(SummarizedExperiment)
library(pSI)

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/pSI_input_Lake_2018_renamed_avg_broadtypes.RData")

#Change NA to 0 for pSI, then make file that filters out anything less than 0.2
avg_lake_2018_all_renamed_broad[is.na(avg_lake_2018_all_renamed_broad)] <- 0
avg_lake_2018_all_renamed_broad_notfiltered <- avg_lake_2018_all_renamed_broad
avg_lake_2018_all_renamed_broad[avg_lake_2018_all_renamed_broad < .2] <- NA

#Need to select more specific p value?
#need to use NA filter? Maybe for lowly expressed or for which genes for each cell type? maybe already done...
#need to do any normalization first?
lake_spi <- specificity.index(avg_lake_2018_all_renamed_broad_notfiltered, avg_lake_2018_all_renamed_broad, p_max=1)

pSI.count(lake_spi)
pSI.list(lake_spi, write.csv=FALSE)		

#load objects

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_turquoise.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_blue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_brown.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_yellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_green.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_red.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_black.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_pink.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_magenta.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_purple.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_greenyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_tan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_salmon.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_cyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_midnightblue.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightcyan.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_grey60.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightgreen.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_lightyellow.rda")
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/rdas/gene_Symbols_WGCNA_royalblue.rda")

#Calculation

fisher.iteration(lake_spi,genes_wgcna_grey)
fisher.iteration(lake_spi,genes_wgcna_turquoise)	
fisher.iteration(lake_spi,genes_wgcna_blue)	
fisher.iteration(lake_spi,genes_wgcna_brown)	
fisher.iteration(lake_spi,genes_wgcna_yellow)	
fisher.iteration(lake_spi,genes_wgcna_green)	
fisher.iteration(lake_spi,genes_wgcna_red)	
fisher.iteration(lake_spi,genes_wgcna_black)	
fisher.iteration(lake_spi,genes_wgcna_pink)	
fisher.iteration(lake_spi,genes_wgcna_magenta)	
fisher.iteration(lake_spi,genes_wgcna_purple)	
fisher.iteration(lake_spi,genes_wgcna_greenyellow)	
fisher.iteration(lake_spi,genes_wgcna_tan)	
fisher.iteration(lake_spi,genes_wgcna_salmon)	
fisher.iteration(lake_spi,genes_wgcna_cyan)	
fisher.iteration(lake_spi,genes_wgcna_midnightblue)	
fisher.iteration(lake_spi,genes_wgcna_lightcyan)	
fisher.iteration(lake_spi,genes_wgcna_grey60)	
fisher.iteration(lake_spi,genes_wgcna_lightgreen)	
fisher.iteration(lake_spi,genes_wgcna_lightyellow)	
fisher.iteration(lake_spi,genes_wgcna_royalblue)	
	
#Save output/results

write.csv(fisher.iteration(lake_spi,genes_wgcna_grey), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey_Lake_broad.csv")				
write.csv(fisher.iteration(lake_spi,genes_wgcna_turquoise), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_turquoise_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_blue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_blue_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_brown), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_brown_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_yellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_yellow_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_green), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_green_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_red), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_red_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_black), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_black_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_pink), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_pink_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_magenta), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_magenta_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_purple), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_purple_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_greenyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_greenyellow_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_tan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_tan_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_salmon), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_salmon_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_cyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_cyan_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_midnightblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_midnightblue_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightcyan), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightcyan_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_grey60), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_grey60_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightgreen), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightgreen_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_lightyellow), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_lightyellow_Lake_broad.csv")		
write.csv(fisher.iteration(lake_spi,genes_wgcna_royalblue), file = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/csvs/Results/genes_WGCNA_royalblue_Lake_broad.csv")		
