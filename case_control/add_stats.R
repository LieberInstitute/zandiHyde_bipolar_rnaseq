###

## load significant hits
load("bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_FDR05either.rda")

## load interaction terms
load("interaction_model_results.rda")
intStats = rbind(outGene_bothRegion[rownames(outGene_bothRegion) %in% rownames(sigOut),3:11],
	outExon_bothRegion[rownames(outExon_bothRegion) %in% rownames(sigOut),3:11],
	outJxn_bothRegion[rownames(outJxn_bothRegion) %in% rownames(sigOut),3:11],
	outTx_bothRegion[rownames(outTx_bothRegion) %in% rownames(sigOut),3:11])
intStats = intStats[rownames(sigOut),]	

