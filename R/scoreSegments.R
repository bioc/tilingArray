##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
scoreSegments = function(x, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  minOverlap = 0.8,
  verbose = TRUE) {

  gff$Name   = getAttributeField(gff$attributes, "Name")
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
        pt                    = rep(as.numeric(NA), cp),
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
        if(verbose)
          cat(wgff, ": ", nrow(sgff), "  ", sep="")
        
        p = function(x) paste(wgff, x, sep=".")

        ov = rep(as.numeric(NA), cp)        
        dl = dr = rep(as.integer(NA), cp)   
        ft = character(cp)  
        
        for(j in 1:cp) {
          ssj = segStart[j]
          sej = segEnd[j]

          ## this is the overlap (between 0 and 1) of every feature with the current segment
          overlap = (pmin(sej, sgff$end) - pmax(ssj, sgff$start) + 1) / sgff$length
          whf = which(overlap > minOverlap)
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
      
      theCut = cut(seq(along=dat$y), th-1, labels=paste("<=", th[-1]-1, sep=""))
      
      ## use only the non-duplicated probes!
      ys     = split(dat$y[dat$xunique], theCut[dat$xunique])
      lls    = listLen(ys)
      means  = sapply(ys, mean)
      sds    = sapply(ys, sd)
      p      = pt(means/sds*sqrt(lls), df=lls-1, lower.tail=FALSE)
      
      ## there might not be data for all cut levels, since some have been
      ## dropped for being duplicated
      mt = match(names(means), levels(theCut))
      stopifnot(identical(names(means), names(sds)),
                identical(names(means), names(p)), !any(is.na(mt)))
      segScore$level[idx[mt]]  = means
      segScore$pt[idx[mt]]     = p
      
      ## just check
      ri = sample(which(listLen(ys) > 5), size=1)
      stopifnot(abs(p[ri] - t.test(ys[[ri]], alternative="greater")$p.value) < 1e-10)
      
      rv = rbind(rv, segScore)
    } ## for strand
  } ## for chr
  return(rv)
}
