##
## categorize segments for the UTR mapping
##
categorizeSegmentsUTRmap = function(env, maxDuplicated=0.5, zThresh=2) {

  s = categorizeSegmentsPie(env, maxDuplicated=maxDuplicated)

  minZ =  pmin(s[,"zLeft"], s[,"zRight"])
  hasGoodFlanks = (!is.na(minZ) & (minZ >= zThresh))

  ## annotated ORF, z-score criterion, nuclear
  sel = ((s[, "category"]=="annotated ORF") & hasGoodFlanks &
         !is.na(s[, "utr3"]) & !is.na(s[, "utr5"]) & (s[,"chr"] <= 16))

  s$goodUTR = ifelse(sel, minZ, as.numeric(NA))

  return(s)
}

##
## Categorize segments for the Pie chart and for the cross-species comparison.
## For each segment, the factor 'category' contains an assignment to a category  
##
categorizeSegmentsPie = function(env, maxDuplicated=0.5,
  minNewSegmentLength=24,
  zThresh=1) {

  ## results data structure: a factor which assigns a category to each segment:
  catg = factor(rep(NA, nrow(s)), 
    levels = c("annotated ORF", "ncRNA", "other annotation",
      "novel isolated", "novel antisense", 
      "untranscribed", "excluded"))

  s = get("segScore", env)
  threshold = get("threshold", env)
  stopifnot(length(threshold)==1, is(s, "data.frame"))

  ## Step 1: frac.dup
  catg[ s[,"frac.dup"] >= maxDuplicated ] = "excluded"

  ## Step 2: untranscribed
  catg[ is.na(catg) & s[,"level"] < threshold ] = "untranscribed"

  ## step 3: annotated
  wh  = which(is.na(catg))
  ovF = strsplit(s[wh,"mostOfFeatureInSegment"],  split=", ")
  
  for(aCateg in levels(catg)[3:1]) {
    theseNames = switch(aCateg,
      ## orf_classification: applies to "gene", "CDS", and "intron"
      "annotated ORF"   = gff$Name[(gff[, "feature"]=="gene") & (gff[, "orf_classification"] %in% c("Verified", "Uncharacterized"))], 
      "ncRNA"           = gff$Name[(gff[, "feature"] %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA"))], 
      "other annotation"= gff$Name[((gff[, "feature"]=="gene") & (gff[, "orf_classification"]=="Dubious")) |
                                    (gff[, "feature"] %in% c("transposable_element", "transposable_element_gene"))],
      stop("Zapperlot")
    )
    stopifnot(!is.null(theseNames))

    ## first with features that this segment is contained in
    isThisCateg = sapply(ovF, function(x) any(x %in% theseNames))
    catg[wh[isThisCateg]] = aCateg

  } ## for aCateg

  ## step 4: overlap>0 but <50% (overlappingFeature)
  catg[ is.na(catg) & (s[, "overlappingFeature"]!="") ] = "other annotation"
  
  ## step 5: novelty filter
  zmin = pmin(s[, "zLeft"], s[, "zRight"])
  catg[ is.na(catg) & ( is.na(zmin) | (zmin <zThresh) |
        (s[,"length"] < minNewSegmentLength) |
        (s[,"oppositeExpression"] > threshold)) ] = "excluded"

  ## step 6: novel - isolated or antisense
  catg[ is.na(catg) & (s[,"oppositeFeature"]=="") ] = "novel isolated"
  catg[ is.na(catg) ] = "novel antisense"

  s$category = catg

  return(s)
}


