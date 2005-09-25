##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
colMedians = function(x)
  apply(x, 2, median)

vectornorm = function(x) {
  sqrt(mean(x*x))
}

zscore = function(x, x0) {
  if(nrow(x)>2) {
    (mean(x0-colMeans(x)))/vectornorm(sd(x))
  } else {
    as.numeric(NA)
  }
}

movingWindow = function(x, y, width) {
  stopifnot(length(x)==nrow(y))
  w = which(x+width-1 <= x[length(x)])
  if(length(w)==0) {
    res = median(y)
  } else {
    res = sapply(w, function(i)
      median(y[x>=x[i] & x<(x[i]+width)]))
    stopifnot(!any(is.na(res)))
    res = min(res)
  }
  res
}

scoreSegments = function(s, gff, 
  nrBasePerSeg = 1500, # average segment length in base-pairs 
  probeLength  = 25, 
  params = c(overlapFraction = 0.5, oppositeWindow = 100, flankProbes=10),
  verbose = TRUE) {

  data(yeastFeatures)
  stopifnot(all(gff$feature %in% rownames(yeastFeatures)))
  transcribedFeatures = rownames(yeastFeatures)[yeastFeatures$isTranscribed]
    
  rv = NULL # initialize result
  for(chr in chrs) { ## chrs (chromosomes) not defined anywhere!
    for(strand in c("+", "-")) {
      
      dat     = get(paste(chr, strand, "dat", sep="."), s)
      datOppo = get(paste(chr, otherStrand(strand), "dat", sep="."), s)
      seg     = get(paste(chr, strand, "seg", sep="."), s)

      stopifnot(is.numeric(dat$unique), is.numeric(datOppo$unique))

      lengthChr = dat[["end"]][length(dat[["end"]])]
      cp        = round(lengthChr/nrBasePerSeg)
      # cp: expected number of segments on this chromosome
      dzz = cp - nrow(seg[["th"]])
      if(dzz>0) {
        if(dzz<=2) {
           cp = nrow(seg[["th"]])
         } else {
           stop("'nrBasePerSeg' is too small for 's'")
         }
      }
      if(verbose)
        cat(chr, ".", strand, ": length=", lengthChr, ", ", cp, " segments. ", sep="")

      segScore = data.frame(
        chr                   = rep(as.integer(NA), cp),
        strand                = I(character(cp)),
        start                 = rep(as.integer(NA), cp),
        end                   = rep(as.integer(NA), cp),
        length                = rep(as.integer(NA), cp),
        level                 = rep(as.numeric(NA), cp),
        featureInSegment      = I(character(cp)),
        mostOfFeatureInSegment= I(character(cp)),
        overlappingFeature    = I(character(cp)),
        overlapFeatAll        = I(character(cp)),
        oppositeFeature       = I(character(cp)),
        oppositeExpression    = rep(as.numeric(NA), cp),
        utr5                  = rep(as.integer(NA), cp),
        utr3                  = rep(as.integer(NA), cp),
        z5                    = rep(as.numeric(NA), cp),
        z3                    = rep(as.numeric(NA), cp),
        dist5                 = rep(as.integer(NA), cp),
        dist3                 = rep(as.integer(NA), cp),
        frac.dup              = rep(as.numeric(NA), cp))
      
      ## th[i] is 1 + (end point of segment i) which is the same as
      ## the start point of segment (i-1).
      th =  c(1, seg[["th"]][cp, 1:cp])
      i1 = th[-length(th)]         ## start points
      i2 = th[-1] - 1              ## end points

      ## extract relevant data from "dat"
      wh = which(dat[["ss"]])
      dStart = dat[["start"]][wh]      ## start base of all probes 
      dEnd   = dat[["end"]][wh]        ## end base of all probes
      dUniq  = dat[["unique"]][wh]
      dY     = dat[["y"]][wh,, drop=FALSE]

      ## extract relevant data from "datOppo"
      wh = which(datOppo[["ss"]])
      dOppoStart = datOppo[["start"]][wh] ## start base of all probes
      dOppoEnd   = datOppo[["end"]][wh]   ## end base of all probes
      dOppoUniq  = datOppo[["unique"]][wh]
      dOppoY     = datOppo[["y"]][wh,, drop=FALSE]
        
      ## double-check: ascending?
      stopifnot(all(diff(dStart+dEnd)>=0), all(diff(dOppoStart+dOppoEnd)>=0))

      ## ... and insert into the results table
      segScore[, "chr"]      = chr
      segScore[, "strand"]   = strand
      segScore[, "start"]    = dStart[i1]
      segScore[, "end"]      = dEnd[i2]
      segScore[, "length"]   = dEnd[i2]-dStart[i1]+1

      sel = (gff[, "chr"]==chr) & (gff[, "feature"] %in% transcribedFeatures)
      same.gff = gff[ sel & gff[, "strand"]==strand, ] 
      oppo.gff = gff[ sel & gff[, "strand"]==otherStrand(strand), ]

      ## for overlapFeatAll: features on either strand
      all.gff = gff[sel, ]

      stopifnot(!any(is.na(all.gff[, "Name"])),
                !any(is.na(same.gff[, "Name"])),
                !any(is.na(oppo.gff[, "Name"])))
                
      ## 'incompleteFeature' is used below, in the definition of "featureInSegment"
      ## Note: for the 050909 segmentation, the following definition was used:
      ## incompleteFeature ="CDS"
      ## this is one is more correct, for subsequent versions:
      incompleteFeature =c("CDS", "intron") ## wh 25.Sep 2005
      stopifnot(incompleteFeature %in% transcribedFeatures)
      
      utrLeft  = utrRight = dl = dr = rep(as.integer(NA), cp)   
      ft1 = ft2 = ft3 = ft4 = ft5 = character(cp)  
      zl = zr = lev = oe = fd = rep(as.numeric(NA), cp)
      nrFlankProbes = params["flankProbes"]
      
      for(j in 1:cp) {
        startj = dStart[i1[j]]
        endj   = dEnd[i2[j]]
        
        ## data from segment, and opposite
        kw     = (dStart>=startj) & (dEnd<=endj)
        ksel   = (dUniq==0) & kw
        ym     = dY[ksel,,drop=FALSE]

        ## frac.dup
        fd[j]  = 1-mean(dUniq[kw]==0)
        
        ksel   = (dOppoUniq==0) & (dOppoStart>=startj) & (dOppoEnd<=endj)
        xOppo  = (dOppoStart[ksel]+dOppoEnd[ksel])/2
        yOppo  = dOppoY[ksel,,drop=FALSE]
        
        cmym   = colMedians(ym)
        lev[j] = median(ym)

        if(is.na(lev[j]))
          stopifnot(fd[j]>0.5)
                  
        ## data from flanks, for segment quality scores
        if(j>1) {
          probesLeft = which((dUniq==0) & (dEnd<startj) & (dStart>=dStart[i1[j-1]]))
          if(length(probesLeft) > nrFlankProbes)
            probesLeft = probesLeft[1:nrFlankProbes]
          yl = dY[probesLeft,, drop=FALSE]
          zl[j] = zscore(yl, cmym)
        }
        if(j<cp) {
          probesRight = which((dUniq==0) & (dStart>endj) & (dEnd<=dEnd[i2[j+1]]))
          if(length(probesRight) > nrFlankProbes)
            probesRight = probesRight[1:nrFlankProbes]
          yr    = dY[probesRight,, drop=FALSE]
          zr[j] = zscore(yr, cmym)
        } 

        ## distance to next features on the left and on the right:
        dl[j]  = posMin(startj - same.gff[, "end"])
        dr[j]  = posMin(same.gff[, "start"] - endj)
        dlOppo = posMin(startj - oppo.gff[, "end"])
        drOppo = posMin(oppo.gff[, "start"] - endj)
        
        nm1 = nm2 = nm3 = nm4 = nm5 = character(0)
        ## featureInSegment: feature is fully contained in the segment
        ## a 'feature' is one of the things in 'known_features', except 'CDS'
        ## (since we want the whole gene, not just its exons)
        wh1 = which(
          (same.gff[, "start"] >= startj) &
          (same.gff[, "end"]   <= endj ) &
          (same.gff[, "feature"] != incompleteFeature))
        if(length(wh1)>0) {
          nm1 = unique(same.gff[wh1, "Name"])
          stopifnot(!any(duplicated(nm1)))
          ft1[j] = paste(nm1, collapse=", ")
        }

        ## mostOfFeatureInSegment: overlap is more than 50% of feature length
        overlapSame   = pmin(endj, same.gff[,"end"]) - pmax(startj, same.gff[,"start"]) + 1 
        wh2 = which( overlapSame / (same.gff[,"end"]-same.gff[,"start"]+1) >=
          params[["overlapFraction"]])
        if(length(wh2)>0) {
          nm2    = unique(same.gff[wh2, "Name"])
          ft2[j] = paste(nm2, collapse=", ")
        }
        stopifnot(all(nm1 %in% nm2))

        ## overlappingFeature: any overlap with the segment
        wh3 = which( overlapSame > 0)
        if(length(wh3)>0) {
          nm3    = unique(same.gff[wh3, "Name"])
          ft3[j] = paste(nm3, collapse=", ")
        }
        stopifnot(all(nm2 %in% nm3))

        ## oppositeFeature 
        overlapOppo = (pmin(endj, oppo.gff[, "end"]) - pmax(startj, oppo.gff[,"start"]))
        wh4         = which( overlapOppo > 0 )
        if(length(wh4)>0)
          ft4[j] = paste(unique(oppo.gff[wh4, "Name"]), collapse=", ")

        ## overlapFeatAll
        wh5 =  which((all.gff[,"start"] <= endj) & (all.gff[, "end"] >= startj))
        if(length(wh5)>0) {
          nm5 = unique(all.gff[wh5, "Name"])
          ft5[j] = paste(nm5, collapse=", ")
        }
        stopifnot(all(nm3 %in% nm5))
        
        ## UTR lengths: if the segment contain exactly one gene
        if((length(wh1)==1) && (same.gff[wh1, "feature"]=="gene") &&
           identical(nm1, nm2) && identical(nm1, nm3)) {
          utrLeft[j]  = same.gff[wh1, "start"] - startj 
          utrRight[j] = endj - same.gff[wh1, "end"]
        }

        ## expression on opposite strand?
        oe[j] = movingWindow(x=xOppo, y=yOppo, width=params[["oppositeWindow"]])

      } ## for j

      stopifnot(!any(is.na(fd)))

      switch(strand,
        "+" = {
          segScore[, "utr5"]  = utrLeft
          segScore[, "utr3"]  = utrRight
          segScore[, "z5"]    = zl
          segScore[, "z3"]    = zr
          segScore[, "dist5"] = dl
          segScore[, "dist3"] = dr
        },
        "-" = {
          segScore[, "utr3"] = utrLeft
          segScore[, "utr5"] = utrRight
          segScore[, "z3"]   = zl
          segScore[, "z5"]   = zr
          segScore[, "dist3"] = dl
          segScore[, "dist5"] = dr
        },
        stop("Zapperlot"))

      segScore[, "featureInSegment"]       = ft1
      segScore[, "mostOfFeatureInSegment"] = ft2
      segScore[, "overlappingFeature"]     = ft3
      segScore[, "oppositeFeature"]        = ft4
      segScore[, "overlapFeatAll"]         = ft5
      segScore[, "oppositeExpression"] = oe
      segScore[, "level"]              = lev
      segScore[, "frac.dup"]           = fd

      ### missing: enter confidence interval for each change-point (right segment border)
        
      rv = rbind(rv, segScore)

      if (verbose)
        cat("\n")
    } ## for strand
  } ## for chromosome
  return(rv)
}
