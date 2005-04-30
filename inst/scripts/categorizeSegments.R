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
  isUnanno = (s$featureInSegment=="" & s$segmentInFeature=="" & isUnique)
  thresh = calcThreshold(s$level, sel=isUnanno, showPlot=FALSE, main=rt)
  cat("thresh=", signif(thresh, 2), "\n")

  isTranscribed = (isUnique & (s$level>=thresh))
  neighborSegmentUnTranscribed = ( c(FALSE, s$level[-nrow(s)]<thresh) &
                                   c(s$level[-1]<thresh, FALSE) )
  
  ll = listLen(strsplit(s$featureInSegment, split=", "))
  isWellDefined = !(is.na(s$sdLeft) | is.na(s$sdRight) | is.na(s$excurse) | ll !=1 )
  isGoodGene    = (s$featureInSegment %in% goodGenes)   ## c

  candidates = which(isTranscribed & neighborSegmentUnTranscribed & isWellDefined & isGoodGene)

  ## max. rank
  maxrk = pmax(rank(s$sdLeft[candidates]), rank(s$sdRight[candidates]),
              rank(s$excurse[candidates]))
     
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

categorizeSegmentsPie = function(s, maxDuplicated=0.5) {
  isUnique = (s$frac.dup < maxDuplicated)
  isUnanno = (s$featureInSegment=="" & s$segmentInFeature=="" & isUnique)
  thresh = calcThreshold(s$level, sel=isUnanno, showPlot=TRUE, main=rt)
  cat("thresh=", signif(thresh, 2), "\n")

  isTranscribed = (isUnique & (s$level>=thresh))
  wh = which(isTranscribed)
  
  category = factor(rep(NA, nrow(s)), 
    levels = c("verified", "ncRNA", "uncharacterized", "dubious", "other",
      "unA", "unI", "not expressed"))
  category[ -wh ] = "not expressed"

  ## split
  patFinS = strsplit(s$featureInSegment[wh],  split=", ")
  patSinF = strsplit(s$segmentInFeature[wh],  split=", ")
  
  count           = matrix(NA, nrow=length(levels(category)), ncol=2)
  rownames(count) = levels(category)
  colnames(count) = c("observed", "in genome")
  count = count[-which(rownames(count)=="not expressed"), ]
  
  for(categ in levels(category)[5:1]) {
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
    count[categ, "observed"]  = length(unique(intersect(theseNames,
           union(unlist(patFinS), unlist(patSinF)))))

    ## first with features that this segment is contained in
    isThisCategSinF = sapply(patSinF, function(x) any(!is.na(match(x, theseNames))))
    category[wh][isThisCategSinF] = categ

    ## then with features contained in this segment
    isThisCategFinS = sapply(patFinS, function(x) any(!is.na(match(x, theseNames))))
    category[wh][isThisCategFinS] = categ

  } ## for categ
  
  ## Properly count the unannotated segments: they should be flanked on both sides by something
  ## which is not expressed. Possible consecutive expressed unannotated segments are merged.
  category[isTranscribed & isUnanno & s$oppositeFeature==""] = "unI"
  category[isTranscribed & isUnanno & s$oppositeFeature!=""] = "unA"
  stopifnot(!any(is.na(category)))

  unBoth = (category %in% c("unA", "unI"))
  n = length(unBoth)
  isConsecutive = c(FALSE, unBoth[2:n] & unBoth[1:(n-1)])
  
  mc  = category[!isConsecutive]

  n   = length(mc)
  s1 = (mc[1:(n-2)]=="not expressed" & mc[2:(n-1)]=="unA" & mc[3:n]=="not expressed")
  s2 = (mc[1:(n-2)]=="not expressed" & mc[2:(n-1)]=="unI" & mc[3:n]=="not expressed")
  count["unA", ] = c(sum(s1), NA)
  count["unI", ] = c(sum(s2), NA)

  cat("New segments: started with ", sum(unBoth), ", merging resulted in ",
      sum(mc %in% c("unA", "unI")), ",\n'neighbor-not-expressed' filter resulted in ",
      count["unA",1], "+", count["unI", 1], "=", count["unA",1]+count["unI", 1], ".\n\n",
      sep="")

  list(category=category, count=count)
}


