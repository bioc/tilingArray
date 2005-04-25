##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
vectornorm = function(x) sqrt(mean(x*x))

isGoodUTRMappingCandidate = function(xleft, x, xright, minN=5) {
  sl = sx = sr = ex = mx = px = as.numeric(NA)
  stopifnot(is.matrix(xleft), is.matrix(x), is.matrix(xright))
  if(nrow(xleft)>=minN)
    sl = vectornorm(sd(xleft[nrow(xleft)-(0:(minN-1)), ]))
  if(nrow(xright)>=minN)
    sr = vectornorm(sd(xright[1:minN, ]))
  if(nrow(x)>=minN) {
    sx = vectornorm(sd(x))
    mx = mean(x)
    k  = 3:nrow(x)
    x  = matrix(colMeans(x), ncol=ncol(x), nrow=nrow(x), byrow=TRUE) - x
    ex = max(rowMeans(x[k-2, ]+x[k-1, ]+x[k, ])) / 3
  }
  return(c(sl, sx, sr, mx, ex))
}

scoreSegments = function(s, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  params = c(minOverlapFractionSame = 0.8, minOverlapBasesOppo = 40),
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
      
      dat  = get(paste(chr, strand, "dat", sep="."), s)
      seg  = get(paste(chr, strand, "seg", sep="."), s)
      cp   = round(max(dat$x)/nrBasePerSeg)
      dzz  = cp - nrow(seg$th)
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
        excurse               = rep(as.numeric(NA), cp),
        sdLeft                = rep(as.numeric(NA), cp),
        sdThis                = rep(as.numeric(NA), cp),
        sdRight               = rep(as.numeric(NA), cp),
        frac.dup              = rep(as.numeric(NA), cp),
        same.feature          = I(character(cp)),
        same.overlap          = rep(as.numeric(NA), cp),
        same.dist5            = rep(as.integer(NA), cp),
        same.dist3            = rep(as.integer(NA), cp),
        oppo.feature          = I(character(cp)),
        oppo.overlap          = rep(as.numeric(NA), cp),
        oppo.dist5            = rep(as.integer(NA), cp),
        oppo.dist3            = rep(as.integer(NA), cp))
      
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
      segScore$length[idx]   = dat$x[i2]-dat$x[i1]+1
      segScore$frac.dup[idx] = mapply(function(h1, h2) {
        z = dat$xunique[h1:h2]
        1-sum(z)/length(z)
      }, i1, i2)
      
      segStart = segScore$start[idx]+probeMiddle
      segEnd   = segScore$end[idx]+probeMiddle

      for(wgff in c("same", "oppo")) {
        gffstrand = switch(wgff,
          same = strand,
          oppo = otherStrand(strand),
          stop("Sapperlot"))
        sgff = gff[ gff$seqname == chrSeqname[chr] &
          gff$strand  == gffstrand &
          gff$feature %in% knownFeatures, ]
        p = function(x) paste(wgff, x, sep=".")

        ov = rep(as.numeric(NA), cp)        
        dl = dr = rep(as.integer(NA), cp)   
        ft = character(cp)  
        
        for(j in 1:cp) {
          ssj = segStart[j]
          sej = segEnd[j]
          overlap = (pmin(sej, sgff$end) - pmax(ssj, sgff$start) + 1) 
          whf = switch(wgff,
            same = {
              which(overlap/sgff$length > params["minOverlapFractionSame"])
            },
            oppo = {
              which(overlap > params["minOverlapBasesOppo"])
            },
            stop("Zapperlot"))
          
          if(length(whf)>0) {
            ## The segment contains one or more features:
            ft[j] = paste(unique(sgff$Name[whf]), collapse=", ")
            ## one number measuring overlap: from start of leftmost feature in whf to end of
            ## rightmost one:
            leftStart = min(sgff$start[whf])
            rightEnd  = max(sgff$end[whf])
            ov[j] = (min(sej, rightEnd)-max(ssj, leftStart)+1) / (rightEnd-leftStart+1)
            dl[j] = leftStart-ssj+1
            dr[j] = sej-rightEnd+1
          } else {
            ## The segment contains no features:
            ## distance to next features on the left and on the right:
            dl[j] = posMin(ssj - sgff$end)
            dr[j] = posMin(sgff$start - sej)
          }
        } ## for j
        segScore[idx, p("overlap")] = ov
        segScore[idx, p(distleft)]  = dl
        segScore[idx, p(distright)] = dr
        segScore[idx, p("feature")] = ft
      } ## for wgff
      if(verbose)
        cat("\n")

      vars = matrix(as.numeric(NA), nrow=5, ncol=cp)
      stopifnot(all(diff(dat$x)>=0))
      for(j in 1:cp) {
        yr = dat$y[dat$xunique & (dat$x >  dat$x[i2[j]]) & (dat$x<=dat$x[i2[j]]+50), , drop=FALSE]
        yl = dat$y[dat$xunique & (dat$x <  dat$x[i1[j]]) & (dat$x>=dat$x[i1[j]]-50), , drop=FALSE]
        ym = dat$y[dat$xunique & (dat$x >= dat$x[i1[j]]) & (dat$x<=dat$x[i2[j]])   , , drop=FALSE]
        vars[, j] = isGoodUTRMappingCandidate(yl, ym, yr)
      }

      segScore$sdLeft[idx]  = vars[1,]
      segScore$sdThis[idx]  = vars[2,]
      segScore$sdRight[idx] = vars[3,]
      segScore$level[idx]   = vars[4,]
      segScore$excurse[idx] = vars[5,]
      rv = rbind(rv, segScore)
    } ## for strand
  } ## for chr
  return(rv)
}
