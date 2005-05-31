##
## categorize segments for the UTR mapping
##
categorizeSegmentsUTRmap = function(env, maxDuplicated=0.5, zThresh=2) {
  s = categorizeSegments(env, maxDuplicated=maxDuplicated)

  minZ =  pmin(s[,"zLeft"], s[,"zRight"])
  hasGoodFlanks = (!is.na(minZ) & (minZ >= zThresh))

  ## annotated ORF, z-score criterion, nuclear
  sel = ((s[, "category"]=="annotated ORF") & hasGoodFlanks &
         !is.na(s[, "utr3"]) & !is.na(s[, "utr5"]) & (s[,"chr"] <= 16))

  s$goodUTR = ifelse(sel, minZ, as.numeric(NA))
  return(s)
}

##
## Categorize segments 
##
categorizeSegments = function(env, maxDuplicated=0.5,
  minNewSegmentLength=24,
  zThresh=1) {

  feat1 = c("transposable_element", "transposable_element_gene", "ncRNA","snoRNA","snRNA", "tRNA", "rRNA")
  feat2 = c("dubious gene", "uncharacterized gene", "verified gene")
  
  ## results data structure: a factor which assigns a category to each segment:
  overlap = factor(rep(NA, nrow(s)),
    levels = c("<50%", ">=50%, <100%", "100%"))

  catg = factor(rep(NA, nrow(s)), 
    levels = c("excluded", "untranscribed",
      feat1, feat2, 
      "novel isolated - filtered", "novel isolated - excluded",
      "novel antisense - filtered", "novel antisense - excluded")) 

  s = get("segScore", env)
  threshold = get("threshold", env)
  stopifnot(length(threshold)==1, is(s, "data.frame"))

  ## Step 1: frac.dup
  sel = s[,"frac.dup"] >= maxDuplicated
  catg[ sel ] = "excluded"
    
  ## Step 2: untranscribed
  sel = (is.na(catg) & s[,"level"] < threshold)
  catg[ sel ] = "untranscribed"
  s$isUnIso = (sel & (s[, "overlappingFeature"]=="") & (s[,"oppositeFeature"]==""))
  
  ## step 3: annotated
  wh  = which(is.na(catg))
  attrName = c("<50%"         = "overlappingFeature",
               ">=50%, <100%" = "mostOfFeatureInSegment",
               "100%"         = "featureInSegment")

#   categIDs = list(
#         "other annotation"= gff$Name[((gff[, "feature"]=="gene") & (gff[, "orf_classification"]=="Dubious")) |
#                                       (gff[, "feature"] %in% c("transposable_element", "transposable_element_gene"))],
#         "ncRNA"           = gff$Name[ (gff[, "feature"] %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA"))], 
#         "annotated ORF"   = gff$Name[ (gff[, "feature"]=="gene") & (gff[, "orf_classification"] %in% c("Verified", "Uncharacterized"))])

  categIDs = vector(mode="list", length = length(feat1)+3)
  names(categIDs) = c(feat1, feat2)
  
  for(f in feat1)
    categIDs[[f]] = gff[ gff[, "feature"]==f, "Name" ]

  sel = gff[, "feature"]=="gene"
  categIDs[["dubious gene"]]         = gff[ sel & gff[, "orf_classification"]=="Dubious", "Name"]
  categIDs[["uncharacterized gene"]] = gff[ sel & gff[, "orf_classification"]=="Uncharacterized", "Name"]
  categIDs[["verified gene"]]        = gff[ sel & gff[, "orf_classification"]=="Verified", "Name"]

  stopifnot(all(listLen(categIDs)>0))

  ## Loop over <50%, 50-100%, 100%:
  for(i in seq(along=attrName)) {
    ovF = strsplit(s[wh, attrName[i]],  split=", ")

    ## Loop over three annotation classes
    for(j in seq(along=categIDs)) {
      ## find features that this segment is contained in
                   sel = sapply(ovF, function(x) any(x %in% categIDs[[j]]))
         catg[wh[sel]] = names(categIDs)[j]
      overlap[wh[sel]] = names(attrName)[i]
    } ## for j
  } ## i

  ## step 4: novelty filter
  zmin = pmin(s[, "zLeft"], s[, "zRight"])
  
  filt1 = (is.na(zmin) | (zmin <zThresh))
  filt2 = (s[,"length"] < minNewSegmentLength)
  filt3 = (s[,"oppositeExpression"] > threshold)
  filt  = (filt1|filt2|filt3)

  ##cat("Step 4 (novelty filter): ", sum(sel), " (zThresh: ",
  ##      sum(sel1), ", length: ", sum(sel2), ", oppositeExpression: ",
  ##      sum(sel3), ")", " -> excluded\n", sep="")

  ## step 5: novel - isolated or antisense
  iso  = (s[,"oppositeFeature"]=="")
  isna = is.na(catg)

  catg[isna &  iso & !filt ] = "novel isolated - filtered"
  catg[isna &  iso &  filt ] = "novel isolated - excluded"
  catg[isna & !iso & !filt ] = "novel antisense - filtered"
  catg[isna & !iso &  filt ] = "novel antisense - excluded"

  stopifnot(!any(is.na(catg)))
  s$category = catg
  s$overlap  = overlap
  return(s)
}


