##
## categorize segments for the UTR mapping
##

## For criterion c
goodGenes = gff$Name[gff$feature=="gene" &
  gff$orf_classification %in% c("Uncharacterized", "Verified")]

splicedGenes1 = sort(gff$Name[gff$feature=="intron"])
splicedGenes2 = names(which(table(gff$Name[gff$feature=="CDS"])>=2))
goodGenes = setdiff(goodGenes, union(splicedGenes1, splicedGenes2))

categorizeSegmentsUTRmap = function(s, minOverlap=1, maxDuplicated=0.5) {
  isUnique = (s$frac.dup < maxDuplicated)        ## a
  isUnanno = (s$same.feature=="")
  thresh   = calcThreshold(s$level, sel=isUnique&isUnanno, main=rt)

  ##  cat(rt, ": thresh=", signif(thresh, 2), "\n", sep="")

  isOneGene  = (listLen(strsplit(s$same.feature, split=", "))==1) & (s$same.overlap >= minOverlap) ## b
  isGoodGene = (s$same.feature %in% goodGenes)   ## c

  k = 2:(nrow(s)-1)  ## d+e
  isTranscribed = rep(FALSE, nrow(s))
  isTranscribed[k] = (s$level[k]>=thresh) & (s$level[k-1]<thresh) & (s$level[k+1]<thresh)

  sel = isUnique & isAnno & isGoodGene & isTranscribed

}

##
## categorize segments for the Pie chart
##

categorizeSegmentsPie = function(name, minOverlap=0.8, maxDuplicated=0.5) {
  select = (s$frac.dup < maxDuplicated)
  unanno = (s$same.feature=="" | s$same.overlap < minOverlap)  & select
  thresh = calcThreshold(s$level, sel=unanno, showPlot=TRUE, main=rt)
    
  ## cat("=========", rt, "==========\n")
  cat("thresh=", signif(thresh, 2), "\n")

  transcribed = (select & (s$level>=thresh))
  name = s$same.feature[transcribed]
  stopifnot(!any(is.na(name)))
  
  care = c("verified", "uncharacterized", "dubious",
           "other ncRNA", "unannotated")

  mt = matrix(0, nrow=length(name), ncol=length(care))
  colnames(mt)=care

  ## convert into regular expression
  pat = strsplit(name,  split=", ")
  pat[listLen(pat)==0] = "Unannotated"
  
  nrInGenome = integer(length(care))
  names(nrInGenome) = care
  for(f in care) {
    theseNames = switch(f,
      ## orf_classification: applies to "gene", "CDS", and "intron"
      "verified"        = gff$Name[gff$feature=="gene" & gff$orf_classification=="Verified"],
      "uncharacterized" = gff$Name[gff$feature=="gene" & gff$orf_classification=="Uncharacterized"],
      "other ncRNA"     = gff$Name[gff$feature %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA")],
      "dubious"         = gff$Name[gff$feature=="gene" & gff$orf_classification=="Dubious"],
      "unannotated"     = "Unannotated",
      stop("Zapperlot")
      )
    stopifnot(!is.null(theseNames))
    nrInGenome[f] = length(unique(theseNames))
    
    mt[, f] = sapply(pat, function(x) 
      !all(is.na(match(x, theseNames))))
  }
  ## clean up mt - set to zero all lower priority entries:
  for(i in 1:(ncol(mt)-1))
    mt[ mt[,i]>0, (i+1):ncol(mt) ] = 0 
    
  stopifnot(all(rowSums(mt)==1))
  list(mt=mt, nrInGenome=nrInGenome)
}


