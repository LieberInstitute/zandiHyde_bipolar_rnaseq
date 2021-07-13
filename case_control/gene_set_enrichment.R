
gene_set_enrichment <- function(gene_list,
                                modeling_results,
                                fdr_cut = 0.05){
  
  tstats <- modeling_results[, grep("[f|t]_", colnames(modeling_results))]
  fdrs <- modeling_results[, grep("adj.P.Val_", colnames(modeling_results))]
  
  enrichTab <-
    do.call(rbind, lapply(seq(along.with = tstats), function(i) {
      layer <- list( up = tstats[, i] > 0 & fdrs[, i] < fdr_cut,
                     down = tstats[, i] < 0 & fdrs[, i] < fdr_cut)
      
      enrichList <- lapply(layer, function(l){
        lapply(gene_list, function(g) {
          tt <-
            table(
              Set = factor(modeling_results$ensembl %in% g, c(FALSE, TRUE)),
              Layer = factor(l, c(FALSE, TRUE))
            )
          fisher.test(tt)
        })
      })
      
      o <- lapply(enrichList, function(el) {
        o <- data.frame(
          OR = vapply(el, "[[", numeric(1), "estimate"),
          Pval = vapply(el, "[[", numeric(1), "p.value"),
          ID = colnames(tstats)[i],
          stringsAsFactors = FALSE)
        o$test <- rownames(o)
        return(o)   
      })
      o <- do.call("rbind", o)
      o$ID <- paste0(o$ID,"_",ss(rownames(o),"\\."))
      return(o)
    }))
  
  rownames(enrichTab) <- NULL
  enrichTab$fdr_cut <- fdr_cut
  enrichTab$ID <- gsub("_"," ",gsub("t_","",enrichTab$ID))
  return(enrichTab)
}
