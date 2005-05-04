##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
vectornorm = function(x) sqrt(mean(x*x))

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
  params = c(minOverlapFractionSame = 0.8, minOverlapOppo = 40,
    oppositeWindow = 100, utrScoreWidth=100),
  verbose = TRUE) {

  ## minOverlapFractionSame: minimal overlap fraction (between 0 and 1) of a feature
  ##   with the current segment

  if(!"Name" %in% names(gff))
    gff$Name   = getAttributeField(gff$attributes, "Name")
  if(!"length" %in% names(gff))
    gff$length = gff$end - gff$start +1

  probeMiddle = (probeLength-1)/2
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
        stop("Sapperlot"))
      
      dat     = get(paste(chr, strand, "dat", sep="."), s)
      datOppo = get(paste(chr, otherStrand(strand), "dat", sep="."), s)
      seg     = get(paste(chr, strand, "seg", sep="."), s)
      cp      = round(max(dat$x)/nrBasePerSeg)
      
      dzz = cp - nrow(seg$th)
      if(dzz>0) {
        if(dzz<=2) {
           cp = nrow(seg$th)
         } else {
           stop("'nrBasePerSeg' is too small for 's'")
         }
      }
      if(verbose)
        cat(chr, ".", strand, ": ", paste(range(dat$x), collapse="..."),
            ", ", cp, " segments. ", sep="")

      segScore = data.frame(
        chr                   = rep(as.integer(NA), cp),
        strand                = I(character(cp)),
        start                 = rep(as.integer(NA), cp),
        end                   = rep(as.integer(NA), cp),
        length                = rep(as.integer(NA), cp),
        level                 = rep(as.numeric(NA), cp),
        geneInSegment         = I(character(cp)),
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
      th =  c(1, seg$th[cp, 1:cp])
      i1 = th[-length(th)]         ## start points
      i2 = th[-1] - 1              ## end points
      
      ## idx: indices of the segments in the result table
      idx = 1:cp  
      segScore$chr[idx]      = chr
      segScore$strand[idx]   = strand
      segScore$start[idx]    = dat$x[i1]
      segScore$end[idx]      = dat$x[i2]
      segScore$length[idx]   = dat$x[i2]-dat$x[i1]
      segScore$frac.dup[idx] = mapply(function(h1, h2) {
        z = dat$xunique[h1:h2]
        1-sum(z)/length(z)
      }, i1, i2)
      
      segStart = segScore$start[idx]+probeMiddle
      segEnd   = segScore$end[idx]+probeMiddle

      same.gff = gff[ gff$seqname == chrSeqname[chr] &
         gff$strand  == strand &
          gff$feature %in% knownFeatures, ]
      
      oppo.gff = gff[ gff$seqname == chrSeqname[chr] &
          gff$strand  == otherStrand(strand) &
          gff$feature %in% knownFeatures, ]
      
      utrLeft  = utrRight = dl = dr = rep(as.integer(NA), cp)   
      ft1 = ft2 = ft3 = character(cp)  
      zl = zr = lev = oe = rep(as.numeric(NA), cp)
      
      stopifnot(all(diff(dat$x)>=0))
        
      for(j in 1:cp) {
        startj = segStart[j]
        endj   = segEnd[j]

        ## data from segment, and opposite
        ym     = dat$y[dat$xunique & (dat$x   >= dat$x[i1[j]]) & (dat$x    <=dat$x[i2[j]]),, drop=FALSE]
        ksel   = datOppo$xunique & (datOppo$x >= dat$x[i1[j]]) & (datOppo$x<=dat$x[i2[j]])
        xOppo  = datOppo$x[ksel]
        yOppo  = datOppo$y[ksel,,drop=FALSE]
        cmym   = colMeans(ym)
        lev[j] = mean(cmym)
        
        ## data from flanks, for segment quality scores
        yr = dat$y[dat$xunique & (dat$x >  dat$x[i2[j]]) &
          (dat$x<=dat$x[i2[j]]+params["utrScoreWidth"]), , drop=FALSE]
        yl = dat$y[dat$xunique & (dat$x <  dat$x[i1[j]]) &
          (dat$x>=dat$x[i1[j]]-params["utrScoreWidth"]), , drop=FALSE]
        zl[j] = zscore(yl, cmym)
        zr[j] = zscore(yr, cmym)

        ## genes that are fully contained in the segment
        whGinS = which( same.gff$feature=="gene" &
          same.gff$start >= startj &
          same.gff$end   <= endj )
        if(length(whGinS)>0) {
          nm1 = unique(same.gff$Name[whGinS])
          stopifnot(!any(duplicated(nm1)))
          ft1[j] = paste(nm1, collapse=", ")
          if(length(whGinS)==1) {
            ## The segment contains exactly one feature:
            utrLeft[j] =  same.gff$start[whGinS] - startj 
            utrRight[j] = -same.gff$end[whGinS]   + endj 
          } 
        } else {
          nm1 = character(0)
        }

        ## features that have overlap with the segment
        overlapSame = ((pmin(endj, same.gff$end) - pmax(startj, same.gff$start)) /
                        pmin(endj-startj, same.gff$end-same.gff$start))
        whSinF = which( overlapSame >= params["minOverlapFractionSame"])
        if(length(whSinF)>0) {
          nm2    = unique(same.gff$Name[whSinF])
          ft2[j] = paste(nm2, collapse=", ")
        }
        stopifnot(all(nm1 %in% nm2))
        
        ## distance to next features on the left and on the right:
        dl[j] = posMin(startj - same.gff$end)
        dr[j] = posMin(same.gff$start - endj)
        
        ## annotated feature on opposite strand?
        overlapOppo = (pmin(endj, oppo.gff$end) - pmax(startj, oppo.gff$start))
        whOppo      = which( overlapOppo > 0 )
        ft3[j]      = paste(unique(oppo.gff$Name[whOppo]), collapse=", ")

        ## expression on opposite strand?
        oe[j] = movingWindow(x=xOppo, y=yOppo, width=params["oppositeWindow"])

      } ## for j

      segScore[idx, c("utr5", "utr3")]  = switch(strand,
                "+" = c(utrLeft, utrRight),
                "-" = c(utrRight, utrLeft))

      segScore$geneInSegment[idx]      = ft1
      segScore$overlappingFeature[idx] = ft2
      segScore$oppositeFeature[idx]    = ft3
      segScore$oppositeExpression[idx] = oe
      segScore$level[idx]              = lev
      segScore$distLeft[idx]           = dl
      segScore$distRight[idx]          = dr
      segScore$zLeft[idx]              = zl
      segScore$zRight[idx]             = zr
      rv = rbind(rv, segScore)

      if(verbose)
        cat("\n")
    } ## for strand
  } ## for chr
  return(rv)
}
