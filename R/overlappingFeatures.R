overlappingFeatures = function(chr, strand, start, end, gff, featureSet, minOverlap=0) {
  n = length(chr)
  stopifnot(is.numeric(chr), is.numeric(start),
            is.numeric(end), length(start)==n, length(end)==n,
            all(start<=end))
  if(missing(strand))
    strand = rep(as.character(NA), n)
  stopifnot(is.character(strand), length(strand)==n)
  
  matchfun = function(chrj, strandj, startj, endj) {
    sel = gff[, "chr"]==chrj
    if(!is.na(strandj))
      sel = sel & (gff[, "strand"]==strandj)
    wh  = which(sel)
    
    overlap = (pmin(endj,gff[wh,"end"])-pmax(startj,gff[wh,"start"])+1)/(gff[,"end"]-gff[,"start"]+1)
    whsel   = (overlap > minOverlap)
    
    if(!missing(featureSet))
      whsel = whsel & (gff[wh, "feature"] %in% featureSet)

    gff[wh[whsel], "Name"]
  }

  mapply(matchfun, chr, strand, start, end)
}
