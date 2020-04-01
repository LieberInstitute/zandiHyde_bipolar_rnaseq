####
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(recount)
library(jaffelab)
library(IRanges)
library(org.Hs.eg.db)
library(clusterProfiler)

## feature maps
load("../data/zandiHypde_bipolar_rseGene_n511.rda")
load("../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../data/zandiHypde_bipolar_rseTx_n511.rda")


## common expression map to annotate
rowRanges(rse_gene)$gencodeID = rownames(rse_gene)
rowRanges(rse_gene)$Class = "InGen"
rowRanges(rse_gene)$Type = "Gene"

rowRanges(rse_exon)$Class = "InGen"
rowRanges(rse_exon)$Type = "Exon"

colnames(mcols(rowRanges(rse_jxn)))[13:14] = c("gencodeID", "Symbol")
rowRanges(rse_jxn)$Type = "Junction"

colnames(mcols(rowRanges(rse_tx)))[c(5,8)] = c("gencodeID", "Symbol")
rowRanges(rse_tx)$Class = "InGen"
rowRanges(rse_tx)$Type = "Transcript"

name = c("gencodeID", "Symbol", "Type", "Class")
exprsMap = rbind(as.data.frame(rowRanges(rse_gene))[,name],
	as.data.frame(rowRanges(rse_exon))[,name],
	as.data.frame(rowRanges(rse_jxn))[,name],
	as.data.frame(rowRanges(rse_tx))[,name])
exprsMap = DataFrame(exprsMap)

## add entrez
exprsMap$EntrezID = rowData(rse_gene)$EntrezID[
	match(exprsMap$gencodeID, rownames(rse_gene))]

################
## load

load("bipolarControl_deStats_byRegion_qSVAjoint.rda", verbose=TRUE)
exprsMap = exprsMap[rownames(statOut),]
statOut = cbind(statOut, exprsMap)
# save(statOut, file = "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda")
# write.csv(statOut, file = "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_all.csv")

## filter to any sig for peter
sigOut = statOut[which(statOut$adj.P.Val_sACC < 0.05 | 
	statOut$adj.P.Val_Amyg < 0.05),]
# save(sigOut, file = "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_FDR05either.rda")
# write.csv(sigOut, file = "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_FDR05either.csv")

# # Amyg
# table(sigOut$Type, sigOut$adj.P.Val_Amyg<0.05)[,2]
      # # Exon       Gene   Junction Transcript
       # # 206        104          7          3
# # sACC
# table(sigOut$Type, sigOut$adj.P.Val_sACC<0.05)[,2]
      # # Exon       Gene   Junction Transcript
       # # 959        666        111         99


### do gene set  ###
amygSig = statOut[statOut$adj.P.Val_Amyg < 0.05,]
amygGenes = split(amygSig$EntrezID, amygSig$Type)
amygGenes = lapply(amygGenes, function(x) unique(as.character(x[!is.na(x)])))
names(amygGenes) = paste0(names(amygGenes), "_Amygdala")
lengths(amygGenes)
      # Exon_Amygdala       Gene_Amygdala   Junction_Amygdala Transcript_Amygdala
                 # 65                  87                   6                   3

saccSig = statOut[statOut$adj.P.Val_sACC < 0.05,]
saccGenes = split(saccSig$EntrezID, saccSig$Type)
saccGenes = lapply(saccGenes, function(x) unique(as.character(x[!is.na(x)])))
names(saccGenes) = paste0(names(saccGenes), "_sACC")
lengths(saccGenes)
      # Exon_sACC       Gene_sACC   Junction_sACC Transcript_sACC
            # 368             592              56              89
			
g = c(amygGenes, saccGenes)

u = statOut$EntrezID
u = as.character(u[!is.na(u)])

go = compareCluster(g, univ = u,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
save(go, file = "go_enrichment_bpdControl_unsigned.rda")

goDf = as.data.frame(go)
goDf = goDf[order(goDf$qvalue),]
goDf$GeneRatio = paste0(" ",goDf$GeneRatio)  ## add space so it doesn't look like a date to Excel
write.csv(goDf, file = "bipolarControl_GOenrich_byRegionType_FDRany.csv")





### do gene set - split directions  ###

statOut$dir_Amyg = ifelse(statOut$t_Amyg>0, "UP","DOWN")
statOut$dir_sACC = ifelse(statOut$t_sACC>0, "UP","DOWN")
statOut$group_Amyg = paste0(statOut$dir_Amyg,"_",statOut$Type)
statOut$group_sACC = paste0(statOut$dir_sACC,"_",statOut$Type)

amygSig = statOut[statOut$adj.P.Val_Amyg<0.05,]
amygGenes = split(amygSig$EntrezID, amygSig$group_Amyg)
amygGenes = lapply(amygGenes, function(x) unique(as.character(x[!is.na(x)])))
names(amygGenes) = paste0(names(amygGenes), "_Amyg")
lengths(amygGenes)
      # DOWN_Exon_Amyg       DOWN_Gene_Amyg   DOWN_Junction_Amyg DOWN_Transcript_Amyg         UP_Exon_Amyg         UP_Gene_Amyg     UP_Junction_Amyg   UP_Transcript_Amyg
                  # 37                   40                    3                    1                   28                   47                    3                    2	

## getting errors, try removing groups with <= 3
# amygGenes = amygGenes[-c(3,4,7,8)]
				  
saccSig = statOut[statOut$adj.P.Val_sACC<0.05,]
saccGenes = split(saccSig$EntrezID, saccSig$group_sACC)
saccGenes = lapply(saccGenes, function(x) unique(as.character(x[!is.na(x)])))
names(saccGenes) = paste0(names(saccGenes), "_sACC")
lengths(saccGenes)
      # DOWN_Exon_sACC       DOWN_Gene_sACC   DOWN_Junction_sACC DOWN_Transcript_sACC         UP_Exon_sACC         UP_Gene_sACC     UP_Junction_sACC   UP_Transcript_sACC
                 # 184                  351                   37                   18                  184                  241                   19                   71

g = c(amygGenes,saccGenes)

u = statOut$EntrezID
u = unique(as.character(u[!is.na(u)]))

go = compareCluster(g, univ = u,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
save(go, file = "go_enrichment_bpdControl_up_down2.rda")

## combine Amyg and sACC
goDf = as.data.frame(go)
goDf = goDf[order(goDf$qvalue),]
goDf$GeneRatio = paste0(" ",goDf$GeneRatio)  ## add space so it doesn't look like a date to Excel
write.csv(goDf, file = "bipolarControl_GOenrich_byRegionType_up_down_FDRany.csv")















lapply(split(goDf, goDf$Cluster),head,3)

#########################
## fetal enrichment #####

