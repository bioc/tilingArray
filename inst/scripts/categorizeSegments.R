##
## categorize segments for the UTR mapping
##
categorizeSegmentsUTRmap = function(env, zThresh=2) {
  s = categorizeSegments(env)

  minZ =  pmin(s[,"zLeft"], s[,"zRight"])
  hasGoodFlanks = (!is.na(minZ) & (minZ >= zThresh))

  ## annotated ORF, z-score criterion, nuclear
  catg = s[, "simpleCatg"]
  stopifnot("annotated ORF" %in% levels(catg))
  sel = ((catg=="annotated ORF") & hasGoodFlanks &
         !is.na(s[, "utr3"]) & !is.na(s[, "utr5"]) & (s[,"chr"] <= 16))

  s$goodUTR = ifelse(sel, minZ, as.numeric(NA))
  return(s)
}

##
## Categorize segments 
##
allncRNA = c("ncRNA","snoRNA","snRNA","tRNA","rRNA")
  
categorizeSegments = function(env, minNewSegmentLength=48, zThresh=1) {

  s = get("segScore", env)
  threshold = get("threshold", env)
  stopifnot(length(threshold)==1, is(s, "data.frame"))

  feat1 = c("verified gene", "uncharacterized gene", "dubious gene")
  feat2 = c("tRNA", "rRNA", "snoRNA","snRNA", "ncRNA","transposable_element_gene", "transposable_element")
  
  ## results data structure: a factor which assigns a category to each segment:
  overlap = factor(rep(NA, nrow(s)),
    levels = c("<50%", ">=50%", "complete"))
  
  catg = factor(rep(NA, nrow(s)), 
    levels = c(feat1, feat2, 
      "novel isolated - filtered", "novel isolated - unassigned",
      "novel antisense - filtered", "novel antisense - unassigned",
      "excluded", "untranscribed"))

  simpleCategories = c("annotated ORF", "ncRNA(all)", 
    "novel isolated - filtered",  "novel isolated - unassigned",
    "novel antisense - filtered", "novel antisense - unassigned",
    "excluded", "untranscribed", "dubious gene")
  
  ## Step 1: frac.dup
  sel = s[,"frac.dup"] >= maxDuplicated
  catg[ sel ] = "excluded"
    
  ## Step 2: untranscribed
  sel = (is.na(catg) & s[,"level"] < threshold)
  catg[ sel ] = "untranscribed"
  s$isUnIso = (sel & (s[, "overlapFeatAll"]==""))
  
  ## step 3: annotated
  wh  = which(is.na(catg))
  attrName = c("<50%"  = "overlappingFeature",
               ">=50%" = "mostOfFeatureInSegment",
               "complete"  = "featureInSegment")

  categIDs = vector(mode="list", length = length(feat1)+length(feat2))
  names(categIDs) = c(feat1, feat2)
  
  for(f in feat2)
    categIDs[[f]] = gff[ gff[, "feature"]==f, "Name" ]

  sel = gff[, "feature"]=="gene"
  categIDs[["dubious gene"]]         = gff[ sel & gff[, "orf_classification"]=="Dubious", "Name"]
  categIDs[["uncharacterized gene"]] = gff[ sel & gff[, "orf_classification"]=="Uncharacterized", "Name"]
  categIDs[["verified gene"]]        = gff[ sel & gff[, "orf_classification"]=="Verified", "Name"]

  stopifnot(all(listLen(categIDs)>0))

  ## Loop over <50%, >=50%, complete:
  for(i in seq(along=attrName)) {
    ovF = strsplit(s[wh, attrName[i]],  split=", ")

    ## Loop over three annotation classes
    for(j in rev(seq(along=categIDs))) {
      ## find features that this segment is contained in
      sel = sapply(ovF, function(x) any(x %in% categIDs[[j]]))
      if(any(sel)) {
        catg[wh[sel]]    = names(categIDs)[j]
        overlap[wh[sel]] = names(attrName)[i]
      }
    } ## for j
  } ## i

  ## step 4: novelty filter
  zmin = pmin(s[, "zLeft"], s[, "zRight"])
  
  filt1 = (is.na(zmin) | (zmin <zThresh))
  filt2 = (s[,"length"] < minNewSegmentLength)
  filt3 = (s[,"oppositeExpression"] > threshold)
  filt4 = (s[,"oppositeExpression"] > pmax(threshold, s[,"level"]-1))
  filt  = (filt1|filt2|filt4)

  ## step 5: novel - isolated or antisense
  iso  = (s[,"oppositeFeature"]=="")
  isna = is.na(catg)

  cat("Novelty filter: Considering ", sum(isna), " segments.\n1. z-scores < ", zThresh, ": ",
    sum(isna&filt1), "\n2. length < ", minNewSegmentLength, ": ",
    sum(isna&filt2), "\n3. oppositeExpression > threshold: ",
    sum(isna&filt3), "\n4. oppositeExpression > max(threshold, segment level - 1): ",
    sum(isna&filt4), "\nRejected by (1 or 2 or 4): ",
    sum(isna&filt), ".\n\n", sep="")

  catg[isna &  iso & !filt ] = "novel isolated - filtered"
  catg[isna &  iso &  filt ] = "novel isolated - unassigned"
  catg[isna & !iso & !filt ] = "novel antisense - filtered"
  catg[isna & !iso &  filt ] = "novel antisense - unassigned"

  stopifnot(!any(is.na(catg)))
  s$category = catg
  s$overlap  = overlap

  ## simpleCategory
  simc = factor(rep(NA, nrow(s)), levels=simpleCategories)
  simc[ catg %in% c("uncharacterized gene", "verified gene")] = "annotated ORF"
  simc[ catg %in% c(allncRNA)]  = "ncRNA(all)"
  for(lev in simpleCategories[-(1:2)])
    simc[ s[,"category"]==lev] = lev
  s$simpleCatg = simc

  ## piechart category
  pieNames = c(A="overlap >=50%", B="overlap <50%",
  C="novel isolated - filtered", D="novel isolated - unassigned",
  E="novel antisense - filtered", F="novel antisense - unassigned")
  
  pc = factor(rep(NA, nrow(s)), levels=pieNames)
  pc[overlap  %in% c(">=50%", "complete")] = "overlap >=50%"
  pc[overlap  %in% c("<50%")] = "overlap <50%"
  stopifnot(all(levels(pc)[3:6] %in% levels(catg)))
  for(k in levels(pc)[3:6])
    pc[catg == k] = k

  s$pieCat = pc
  
  return(s)
}


