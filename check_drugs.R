###
library(LIBDpheno)
library(SummarizedExperiment)

## load data
load("data/zandiHypde_bipolar_rseGene_n511.rda")

## get brnums
pheno = toxicant[[1]]
pheno = pheno[match(unique(rse_gene$BrNum), paste0("Br", pheno$brnumerical)),]

## just bipolar
pheno = pheno[which(pheno$primarydx == "Bipolar"),]

## checks
drugs = pheno[,c("brnumerical","lithium", "lifetime_lithium", "antipsychotics",
	"lifetime_antipsych", "antidepressants_ssris","lifetime_antidepress",
	"anticonvulsants","lifetime_anticonvulsant","last_cpze",
	"lifetime_medications_list")]
rownames(drugs) = paste0("Br", drugs$brnumerical)
drugs$brnumerical = NULL

drugs$valproic = grepl("(valp)|(Depakote)",drugs$lifetime_medications_list)
table(drugs$valproic)

## two way
table(drugs$valproic, drugs$lifetime_lithium, 
	dnn = c("Life_Valproic", "Life_Lithium"),useNA="ifany")

drugs$lifetime_medications_list[drugs$valproic]
drugs$BrNum = rownames(drugs)
write.table(drugs, file = "drug_info.tsv", sep="\t",row.names=FALSE)