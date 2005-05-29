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

  audit=TRUE
  
  ## Step 1: frac.dup
  sel = s[,"frac.dup"] >= maxDuplicated
  catg[ sel ] = "excluded"
  if(audit)
    cat("Step 1 (frac.dup):", sum(sel), "-> excluded\n")
    
  ## Step 2: untranscribed
  sel = (is.na(catg) & s[,"level"] < threshold)
  catg[ sel ] = "untranscribed"
  s$isUnIso = (sel & (s[, "overlappingFeature"]=="") & (s[,"oppositeFeature"]==""))
  if(audit)
    cat("Step 2 (level):   ", sum(sel), "-> untranscribed\n")
  
  ## step 3: annotated
  wh  = which(is.na(catg))
  ovF = strsplit(s[wh,"mostOfFeatureInSegment"],  split=", ")
  
  for(aCateg in levels(catg)[3:1]) {
    theseNames = switch(aCateg,
      ## orf_classification: applies to "gene", "CDS", and "intron"
      "annotated ORF"   = gff$Name[ (gff[, "feature"]=="gene") & (gff[, "orf_classification"] %in% c("Verified", "Uncharacterized"))], 
      "ncRNA"           = gff$Name[ (gff[, "feature"] %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA"))], 
      "other annotation"= gff$Name[((gff[, "feature"]=="gene") & (gff[, "orf_classification"]=="Dubious")) |
                                    (gff[, "feature"] %in% c("transposable_element", "transposable_element_gene"))],
      stop("Zapperlot")
    )
    stopifnot(!is.null(theseNames))

    ## first with features that this segment is contained in
    isThisCateg = sapply(ovF, function(x) any(x %in% theseNames))
    catg[wh[isThisCateg]] = aCateg

    if(audit)
      cat("Step 3 (annotated): ", sum(isThisCateg), " -> ", aCateg, "\n", sep="")

  } ## for aCateg

  ## step 4: overlap>0 but <50% (overlappingFeature)
  ## We need this criterion, other we will get many cases
  ## where a feature has been chopped up into several segments
  sel = is.na(catg) & (s[, "overlappingFeature"]!="")
  catg[ sel ] = "other annotation"
  if(audit)
    cat("Step 4 (overlappingFeature):", sum(sel), "-> other annotation\n")
  
  ## step 5: novelty filter
  zmin = pmin(s[, "zLeft"], s[, "zRight"])
  
  sel1 = (is.na(catg) & ( is.na(zmin) | (zmin <zThresh)))
  ## sel1 = (is.na(catg) & ( is.na(zmin) | (zmin <0)))
  sel2 = (is.na(catg) & (s[,"length"] < minNewSegmentLength))
  ## sel3 = (is.na(catg) & (s[,"oppositeExpression"] > threshold))
  sel3 = (is.na(catg) & (s[,"oppositeExpression"] > Inf))
  sel  = sel1|sel2|sel3
  catg[ sel ] = "excluded"
  if(audit)
    cat("Step 5 (novelty filter): ", sum(sel), " (zThresh: ",
        sum(sel1), ", length: ", sum(sel2), ", oppositeExpression: ",
        sum(sel3), ")", " -> excluded\n", sep="")

  ## step 6: novel - isolated or antisense
  sel = is.na(catg) & (s[,"oppositeFeature"]=="")
  catg[ sel ] = "novel isolated"
  if(audit)
    cat("Step 6:", sum(sel), "-> novel isolated\n")

  sel = is.na(catg)
  catg[ sel ] = "novel antisense"
  if(audit)
    cat("Step 6:", sum(sel), "-> novel antisense\n")

  s$category = catg

  return(s)
}


