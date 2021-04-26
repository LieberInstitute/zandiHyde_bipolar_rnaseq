
leveneFast <- function(dat, groups) {
  require(matrixStats)
  N = length(groups)
  k = length(unique(groups))
  
  gIndexes = split(seq(along=groups),groups)
  gLens = as.numeric(sapply(gIndexes,length))
  
  zij = dat
  for(i in seq(along=gIndexes)) {
    ind = gIndexes[[i]]
    zij[,ind] = abs(dat[,ind] - rowMedians(dat[,ind]))
  }
  
  zdotdot = rowMeans(zij)
  zidot = sapply(gIndexes, function(x) rowMeans(zij[,x]))
  
  num = rowSums(gLens*(zidot-zdotdot)^2)
  
  out = zidot
  for(i in seq(along=gIndexes)) {
    ind = gIndexes[[i]]
    out[,i] = rowSums((zij[,ind] - zidot[,i])^2)
  }
  denom = rowSums(out)
  
  W = ((N-k)/(k-1))*(num/denom)
  pval = pf(W, df1 = k-1, df2 = N-k, lower.tail = FALSE)
  
  return(data.frame(cbind(W,pval)))
}
