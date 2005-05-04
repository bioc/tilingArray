##
## categorize segments for the UTR mapping
##

## For criterion c
goodGenes = gff$Name[gff$feature=="gene" &
  gff$orf_classification %in% c("Uncharacterized", "Verified")]

splicedGenes1 = sort(gff$Name[gff$feature=="intron"])
splicedGenes2 = names(which(table(gff$Name[gff$feature=="CDS"])>=2))
goodGenes     = setdiff(goodGenes, union(splicedGenes1, splicedGenes2))

categorizeSegmentsUTRmap = function(s, maxDuplicated=0.5) {
  isUnique = (s$frac.dup < maxDuplicated)
  isUnanno = (s$overlappingFeature=="" & isUnique)
  thresh = calcThreshold(s$level, sel=isUnanno, showPlot=FALSE, main=rt)
  cat("thresh=", signif(thresh, 2), "\n")

  isTranscribed = (isUnique & (s$level>=thresh))
  neighborSegmentUnTranscribed = ( c(FALSE, s$level[-nrow(s)]<thresh) &
                                   c(s$level[-1]<thresh, FALSE) )
  
  ll = listLen(strsplit(s$geneInSegment, split=", "))
  isWellDefined = !(is.na(s$sdLeft) | is.na(s$sdRight) | is.na(s$excurse) | ll !=1 )
  isGoodGene    = (s$geneInSegment %in% goodGenes)   ## c

  candidates = which(isTranscribed & neighborSegmentUnTranscribed & isWellDefined & isGoodGene)

  ## max. rank
  maxrk = pmax(rank(s$sdLeft[candidates]), rank(s$sdRight[candidates]),
              rank(s$excurse[candidates]))
     browser()
  res = rep(as.integer(NA), nrow(s))
  res[candidates] = maxrk
  return(res)
}

##
## categorize segments for the Pie chart and for the cross-species comparison
##
## The factor 'category' contains for each segment an assignment.
## The matrix 'count' contains the number of UNIQUE occurences. For genes and
## ncRNA, uniqueness is defined by name. For unannotated segments, it is 
## defined by non-consecutiveness.

categorizeSegmentsPie = function(s, maxDuplicated=0.5, minNewSegmentLength=24,
  zTresh=1) {

  isUnique = (s$frac.dup < maxDuplicated)
  isUnanno = (s$overlappingFeature=="" & isUnique)
  threshold = calcThreshold(s$level, sel=isUnanno, showPlot=TRUE, main=rt)
  cat("threshold=", signif(threshold, 2), "\n")

  isTranscribed = (isUnique & (s$level>=threshold))
  wh = which(isTranscribed)

  ## results data structure: a factor which assigns a category to each segment:
  s$category = factor(rep(NA, nrow(s)), 
    levels = c("verified", "ncRNA", "uncharacterized", "dubious", "other",
      "unIso", "unAnti", "unDubious", "not expressed"))

  count           = matrix(NA, nrow=length(levels(s$category)), ncol=2)
  rownames(count) = levels(s$category)
  colnames(count) = c("observed", "in genome")
  
  ## >>> Phase 1: everything which is not expressed
  s$category[ -wh ] = "not expressed"

  ## >>> Phase 2: segments that overlap annotated features
  ovF = strsplit(s$overlappingFeature[wh],  split=", ")
  
  for(categ in levels(s$category)[5:1]) {
    theseNames = switch(categ,
      ## orf_classification: applies to "gene", "CDS", and "intron"
      "verified"        = gff$Name[gff$feature=="gene" & gff$orf_classification=="Verified"],
      "uncharacterized" = gff$Name[gff$feature=="gene" & gff$orf_classification=="Uncharacterized"],
      "ncRNA"           = gff$Name[gff$feature %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA")],
      "dubious"         = gff$Name[gff$feature=="gene" & gff$orf_classification=="Dubious"],
      "other"           = gff$Name[gff$feature %in% c("pseudogene", "transposable_element",
                              "transposable_element_gene")],
      stop("Zapperlot")
    )
    stopifnot(!is.null(theseNames))
    count[categ, "in genome"] = length(unique(theseNames))

    ## first with features that this segment is contained in
    isThisCateg = sapply(ovF, function(x) any(x %in% theseNames))
    s$category[wh][isThisCateg] = categ

  } ## for categ

  ## >>> Phase 3: remaining segments 
  isUnassigned = is.na(s$category)
  stopifnot(all((isTranscribed & isUnanno) == isUnassigned))  ## just check
  nua1 = sum(isUnassigned) ## audit
    
  ## 3a: merge adjacent ones:
  ## ... prepare
  n            = length(isUnassigned)
  diffChrLeft  = c(TRUE, s$chr[-1]!=s$chr[-n])
  diffChrRight = c(s$chr[-1]!=s$chr[-n], TRUE)
  unStart      = which(isUnassigned & (c(TRUE, !isUnassigned[-n]) | diffChrLeft))
  unEnd        = which(isUnassigned & (c(!isUnassigned[-1], TRUE) | diffChrRight))
  stopifnot(length(unStart)==length(unEnd), all(unEnd>=unStart))

  ## ... merge
  keep = rep(TRUE, n)
  for(j in which(unEnd>unStart)) {
    i1 = unStart[j]
    i2 = unEnd[j]
    s$level[i1] = sum(s$level[i1:i2]*s$length[i1:i2])/sum(s$length[i1:i2])
    keep[ (i1+1) : i2 ] = FALSE
  }       
  s$end[unStart]       = s$end[unEnd]
  s$length[unStart]    = s$end[unStart]-s$start[unStart]
  s$zRight[unStart]    = s$zRight[unEnd]
  s$distRight[unStart] = s$distRight[unEnd]

  s = s[keep, ]
  nua2 = sum(is.na(s$category))
  
  ## 3b. Require large z-scores on both sides
  ## 3c. Length requirement: 3 probes (24 bases)

  cat("\n\n\n\nATTENTION z is the wrong way round !!!! \n\n")
  s$category[ is.na(s$category) & ((s$length<minNewSegmentLength) | (s$zLeft > -zTresh) |
             (s$zRight > -zTresh)) ] = "other"
##  s$category[ is.na(s$category) & ((s$length<minNewSegmentLength) | (s$zLeft < zTresh) |
##             (s$zRight < zTresh)) ] = "other"
             
  ## >>> Phase 4: assign to "unIso", "unAnti", or "unDubious"
  s$category[ is.na(s$category) & s$oppositeFeature=="" ] = "unIso"
  s$category[ is.na(s$category) & s$oppositeExpression >= threshold] = "unDubious"
  s$category[ is.na(s$category) ] = "unAnti"

  tab = table(s$category)
  count[names(tab), "observed"]  = tab
  count = count[-which(rownames(count) %in% c("other", "not expressed")), ]

  nua3 = count["unAnti",1]
  nua4 = count["unIso",1]
  cat("New segments: started with ", nua1, ", merging resulted in ",
      nua2, ",\n'neighbor-not-expressed' and 'length>=", minNewSegmentLength,
      "' filter resulted in ", nua3, "+", nua4, "=", nua3+nua4, ".\n\n",
      sep="")
  
  list(s=s, count=count, threshold=threshold)
}


