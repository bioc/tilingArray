##
## categorize segments for the UTR mapping
##

## For criterion c
goodGenes = gff$Name[gff$feature=="gene" &
  gff$orf_classification %in% c("Uncharacterized", "Verified")]

splicedGenes1 = sort(gff$Name[gff$feature=="intron"])
splicedGenes2 = names(which(table(gff$Name[gff$feature=="CDS"])>=2))
goodGenes     = setdiff(goodGenes, union(splicedGenes1, splicedGenes2))

categorizeSegmentsUTRmap = function(env, maxDuplicated=0.5, zThresh=2) {
  s = get("segScore", env)
  threshold = get("threshold", env)
  stopifnot(length(threshold)==1)

  minZ =  pmin(s[,"zLeft"], s[,"zRight"])

  isTranscribed = ((s[,"frac.dup"] < maxDuplicated) & (s[,"level"]>=threshold))
  hasGoodFlanks = (minZ >= zThresh)
  isWellDefined = (listLen(strsplit(s[,"geneInSegment"], split=", "))==1)
  isGoodGene    = (s[,"geneInSegment"] %in% goodGenes)   ## c

  candidates = which(isTranscribed & hasGoodFlanks & isWellDefined & isGoodGene)

  s$goodUTR = rep(as.integer(NA), nrow(s))
  s[candidates, "goodUTR"] = minZ[candidates]

  return(s)
}

##
## Categorize segments for the Pie chart and for the cross-species comparison.
## For each segment, the factor 'category' contains an assignment to a category  

categorizeSegmentsPie = function(env, maxDuplicated=0.5,
  minNewSegmentLength=24, zThresh=1) {

  s = get("segScore", env)
  threshold = get("threshold", env)
  stopifnot(length(threshold)==1, is(s, "data.frame"))

  isTranscribed = ((s[,"frac.dup"] < maxDuplicated) & (s[,"level"]>=threshold))
  wh = which(isTranscribed)

  ## results data structure: a factor which assigns a category to each segment:
  catg = factor(rep(NA, nrow(s)), 
    levels = c("verified", "ncRNA", "uncharacterized", "dubious", "other",
      "unIso", "unAnti", "unAnti-Dubious", "not expressed"))

  count           = matrix(NA, nrow=length(levels(catg)), ncol=2)
  rownames(count) = levels(catg)
  colnames(count) = c("observed", "in genome")
  
  ## >>> Phase 1: everything which is not expressed
  catg[-wh] = "not expressed"

  ## >>> Phase 2: segments that overlap annotated features
  ovF = strsplit(s[wh,"overlappingFeature"],  split=", ")
  
  for(aCateg in levels(catg)[5:1]) {
    theseNames = switch(aCateg,
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
    count[aCateg, "in genome"] = length(unique(theseNames))

    ## first with features that this segment is contained in
    isThisCateg = sapply(ovF, function(x) any(x %in% theseNames))
    catg[wh[isThisCateg]] = aCateg

  } ## for aCateg

  ## >>> Phase 3: remaining segments 
  isUnassigned = is.na(catg)
  nua1 = sum(isUnassigned) ## audit
    
  ## 3a: merge adjacent ones:
  ## ... prepare
  n            = length(isUnassigned)
  diffChrLeft  = c(TRUE, s[-1, "chr"] != s[-n, "chr"])
  diffChrRight = c(s[-1, "chr"] != s[-n, "chr"], TRUE)
  unStart      = which(isUnassigned & (c(TRUE, !isUnassigned[-n]) | diffChrLeft))
  unEnd        = which(isUnassigned & (c(!isUnassigned[-1], TRUE) | diffChrRight))
  stopifnot(length(unStart)==length(unEnd), all(unEnd>=unStart))

  ## ... merge
  drop = rep(FALSE, n)
  for(j in which(unEnd>unStart)) {
    i1 = unStart[j]
    i2 = unEnd[j]
    s[i1, "level"] = sum(s[i1:i2, "level"] * s[i1:i2, "length"]) / sum(s[i1:i2, "length"])
    drop[ (i1+1) : i2 ] = TRUE
    s[i1, "isIsolatedSame"] = all(s[i1:i2, "isIsolatedSame"])
    s[i1, "isIsolatedOppo"] = all(s[i1:i2, "isIsolatedOppo"])
  }       
  s[unStart, "end"]       = s[unEnd,   "end"]
  s[unStart, "length"]    = s[unStart, "end"] - s[unStart, "start"]
  s[unStart, "zRight"]    = s[unEnd,   "zRight"]
  s[unStart, "distRight"] = s[unEnd,   "distRight"]

  catg[drop] = "other"
  nua2 = sum(is.na(catg))
  
  ## 3b. Length requirement: 3 probes (24 bases)
  catg[is.na(catg) & (s[,"length"]<minNewSegmentLength)] = "other"
  nua3 = sum(is.na(catg))

  ## 3c. Require large z-scores on both sides
  catg[ is.na(catg) & ((s[,"zLeft"] < zThresh)|(s[,"zRight"] < zThresh)) ] = "other"
             
  ## >>> Phase 4: assign to "unIso", "unAnti", or "unAnti-dubious"
  stopifnot(!any(s[,"isIsolatedOppo"] & s[,"oppositeFeature"]!=""))
  catg[ is.na(catg) & s[,"isIsolatedSame"] & s[,"isIsolatedOppo" ] ] = "unIso"
  catg[ is.na(catg) & s[,"isIsolatedSame"] & (s[,"oppositeFeature"]!="") & (s[,"oppositeExpression"] < threshold)] = "unAnti"
  catg[ is.na(catg) & s[,"isIsolatedSame"] & (s[,"oppositeFeature"]!="") ] = "unAnti-dubious"
  catg[ is.na(catg)  ] = "other"

  tab = table(catg)
  count[names(tab), "observed"]  = tab
  count = count[-which(rownames(count) %in% c("other", "not expressed")), ]

  nua4 = count["unAnti",1]
  nua5 = count["unIso",1]
  cat("New segments: started with ", nua1, ", merging resulted in ",
      nua2, ",\n'length>=", minNewSegmentLength, "' filter resulted in ", nua3,
      ",\n'zLeft, zRight>=", zThresh, "' filter resulted in ",
      nua4, "+", nua5, "=", nua4+nua5, ".\n\n",
      sep="")
  
  s$category = catg
  list(s=s, count=count)
}


