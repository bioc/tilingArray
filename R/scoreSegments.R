##------------------------------------------------------------
## Copyright (2005) Wolfgang Huber
##------------------------------------------------------------
scoreSegments = function(x, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  verbose = TRUE) {

  probeMiddle = (probeLength-1)/2
  rv = NULL
  for(chr in chrs) {
    for(strand in c("+", "-")) {
      dat  = get(paste(chr, strand, "dat", sep="."), s)
      seg  = get(paste(chr, strand, "seg", sep="."), s)
      cp   = round(max(dat$x)/nrBasePerSeg)
      if(verbose)
        cat(chr, ".", strand, ": ", paste(range(dat$x), collapse="..."),
            ", ", cp, " segments. ", sep="")

      segScore = data.frame(
        chr                   = integer(cp),
        strand                = I(character(cp)),
        start                 = integer(cp),
        end                   = integer(cp),
        length                = integer(cp),
        level                 = numeric(cp),
        pt                    = numeric(cp),
        frac.dup              = numeric(cp),
        same.feature          = I(character(cp)),
        same.overlap          = numeric(cp),
        same.dist.start2feat  = integer(cp),
        same.dist.end2feat    = integer(cp),
        oppo.feature          = I(character(cp)),
        oppo.overlap          = numeric(cp),
        oppo.dist.start2feat  = integer(cp),
        oppo.dist.end2feat    = integer(cp))
      
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
      
      for(wgff in c("same", "oppo")) {
        gffstrand = switch(wgff,
          same = strand,
          oppo = otherStrand(strand),
          stop("Sapperlot"))
        sgff = gff[ gff$seqname == chrSeqname[chr] &
          gff$strand  == gffstrand &
          gff$feature %in% knownFeatures, ]
        sgff$Name = getAttributeField(sgff$attributes, "Name")
        if(verbose)
          cat(wgff, ": ", nrow(sgff), "  ", sep="")
        
        ## matchProbes2Feats is a matrix of probes (dat$x) times features
        ## and a matrix element [p,j] is TRUE iff probe p is part of
        ## feature j. A probe is considered part of a feature if its 13th
        ## nucleotide falls into it. That's an arbitrary rule, which is
        ## intended to err on the side of assigning probes to features.
        matchProbes2Feats = isWithinInterval(dat$x+probeMiddle, sgff$start, sgff$end)
        
        p = function(x) paste(wgff, x, sep=".")
        for(j in 1:cp) {
          ## this matrix has as many rows are there are probes in the segment,
          ## and as many columns as there as features
          matchSeg2Feats = matchProbes2Feats[i1[j]:i2[j], ]
          ## features that have overlap with this segment:
          whf = which(colSums(matchSeg2Feats) > 0)
          segScore[j, p("feature")] = paste(unique(sgff$Name[whf]), collapse=", ")
          ## fraction of probes that have overlap with any feature:
          segScore[j, p("overlap")] = mean(rowSums(matchSeg2Feats) > 0)
          ## see man page!
          sp  = segScore$start[j] + probeMiddle
          whs = which(matchSeg2Feats[1, ])
          if(length(whs)>0) {
            segScore[j, p("dist.start2feat")] = min(sp - sgff$start[whs])
          } else {
            segScore[j, p("dist.start2feat")] = posMin(sp - sgff$end)
          }
          ## see man page!
          sp  = segScore$end[j] + probeMiddle
          whe = which(matchSeg2Feats[nrow(matchSeg2Feats), ])
          if(length(wh)>0) {
            segScore[j, p("dist.end2feat")] = min(sgff$end[whe] - sp)
          } else {
            segScore[j, p("dist.end2feat")] = posMin(sgff$start - sp)
          }
        } ## for j
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
      
      ## there might not be data for all cut levels, since some have been dropped for being
      ## duplicated
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
