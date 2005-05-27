##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
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

movingWindow =function(x, y, width) {
  stopifnot(length(x)==nrow(y))
  w = which(x+width-1 <= x[length(x)])
  if(length(w)==0) {
    res = +Inf
  } else {
    res = sapply(w, function(i) {
      rg = which(x>=x[i] & x<(x[i]+width))
      mean(y[rg])
    })
    stopifnot(!any(is.na(res)))
    ## plot(x[w], res); browser()
    res = min(res)
  }
  res
}

scoreSegments = function(s, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA", 
        "transposable_element", "transposable_element_gene"),
  params = c(overlapFraction = 0.5, oppositeWindow = 100, utrScoreWidth=100),
  verbose = TRUE) {

  rv = NULL
  for(chr in chrs) {
    for(strand in c("+", "-")) {
      switch(strand,
        "+" = {
          distleft="dist5"; distright="dist3"
        },
        "-" = {
          distleft="dist3"; distright="dist5"
        },
        stop("Sapperlot")
      ) ## end of switch
      
      dat     = get(paste(chr, strand, "dat", sep="."), s)
      datOppo = get(paste(chr, otherStrand(strand), "dat", sep="."), s)
      seg     = get(paste(chr, strand, "seg", sep="."), s)

      lengthChr = dat[["end"]][length(dat[["end"]])]
      cp        = round(lengthChr/nrBasePerSeg)
      
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
        oppositeFeature       = I(character(cp)),
        oppositeExpression    = rep(as.numeric(NA), cp),
        utr5                  = rep(as.integer(NA), cp),
        utr3                  = rep(as.integer(NA), cp),
        distLeft              = rep(as.integer(NA), cp),
        distRight             = rep(as.integer(NA), cp),
        zLeft                 = rep(as.numeric(NA), cp),
        zRight                = rep(as.numeric(NA), cp),
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
      segScore[, "frac.dup"] = mapply(function(h1, h2) {
        1 - mean(dUniq[h1:h2])
      }, i1, i2)

      same.gff = gff[ gff[, "chr"]==chr & gff[, "strand"]==strand &
          gff[, "feature"] %in% knownFeatures, ]
      
      oppo.gff = gff[ gff[, "chr"] == chr & gff[, "strand"]==otherStrand(strand) &
          gff[, "feature"] %in% knownFeatures, ]
      
      utrLeft  = utrRight = dl = dr = rep(as.integer(NA), cp)   
      ft1 = ft2 = ft3 = ft4 = character(cp)  
      zl = zr = lev = oe = rep(as.numeric(NA), cp)
      
      for(j in 1:cp) {
        startj = dStart[i1[j]]
        endj   = dEnd[i2[j]]

        ## data from segment, and opposite
        ksel   = dUniq & (dStart>=startj) & (dEnd<=endj)
        ym     = dY[ksel,,drop=FALSE]
        
        ksel   = dOppoUniq & (dOppoStart>=startj) & (dOppoEnd<=endj)
        xOppo  = (dOppoStart[ksel]+dOppoEnd[ksel])/2
        yOppo  = dOppoY[ksel,,drop=FALSE]
        
        cmym   = colMeans(ym)
        lev[j] = mean(cmym)
        
        ## data from flanks, for segment quality scores
        Ll = Lr = params[["utrScoreWidth"]]
        if(j>1) {
          Ll = min(Ll, segScore[j-1, "length"])
        }
        if(j<cp) {
          Lr = min(Lr, segScore[j+1, "length"])
        }

        yr = dY[ dUniq & (dStart>endj)   & (dEnd  <=endj+Lr),, drop=FALSE]
        yl = dY[ dUniq & (dEnd  <startj) & (dStart>=startj-Ll),, drop=FALSE]
        zl[j] = zscore(yl, cmym)
        zr[j] = zscore(yr, cmym)

        ## distance to next features on the left and on the right:
        dl[j]  = posMin(startj - same.gff[, "end"])
        dr[j]  = posMin(same.gff[, "start"] - endj)
        dlOppo = posMin(startj - oppo.gff[, "end"])
        drOppo = posMin(oppo.gff[, "start"] - endj)
        
        nm1 = nm2 = nm3 = nm4 = character(0)
        ## featureInSegment: fully contained in the segment
        whFinS = which(
          (same.gff[, "start"] >= startj) &
          (same.gff[, "end"]   <= endj ))
        if(length(whFinS)>0) {
          nm1 = unique(same.gff[whFinS, "Name"])
          stopifnot(!any(duplicated(nm1)))
          ft1[j] = paste(nm1, collapse=", ")
          if(length(whFinS)==1) {
            if(same.gff[whFinS, "feature"]=="gene") {
              ## The segment contains exactly one gene
              utrLeft[j]  = same.gff[whFinS, "start"] - startj 
              utrRight[j] = endj - same.gff[whFinS, "end"]
            }
          } 
        }

        ## mostOfFeatureInSegment: more than 50% (="overlapFraction") overlap with the segment
        overlapSame   = pmin(endj, same.gff[,"end"]) - pmax(startj, same.gff[,"start"]) + 1 
        wh2 = which( overlapSame / (same.gff[,"end"]-same.gff[,"start"]+1) >= params[["overlapFraction"]])
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

        ## expression on opposite strand?
        oe[j] = movingWindow(x=xOppo, y=yOppo, width=params[["oppositeWindow"]])

      } ## for j

      segScore[, c("utr5", "utr3")]  = switch(strand,
                "+" = c(utrLeft, utrRight),
                "-" = c(utrRight, utrLeft))
      segScore[, "featureInSegment"]       = ft1
      segScore[, "mostOfFeatureInSegment"] = ft2
      segScore[, "overlappingFeature"]     = ft3
      segScore[, "oppositeFeature"]        = ft4
      segScore[, "oppositeExpression"] = oe
      segScore[, "level"]              = lev
      segScore[, "distLeft"]           = dl
      segScore[, "distRight"]          = dr
      segScore[, "zLeft"]              = zl
      segScore[, "zRight"]             = zr
        
      rv = rbind(rv, segScore)

      if(verbose)
        cat("\n")
    } ## for strand
  } ## for chr
  return(rv)
}
