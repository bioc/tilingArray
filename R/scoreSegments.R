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
        level                 = rep(as.numeric(NA), cp),
        pt                    = rep(as.numeric(NA), cp),
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
      segScore$length[idx]   = dat$x[i2]-dat$x[i1]+1
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
        sgff$length = sgff$end - sgff$start +1
        if(verbose)
          cat(wgff, ": ", nrow(sgff), "  ", sep="")

        
        p = function(x) paste(wgff, x, sep=".")
        
        for(j in 1:cp) {
          segStart = segScore$start[idx[j]]+probeMiddle
          segEnd   = segScore$end[idx[j]]+probeMiddle

          overlap = (pmin(segEnd, sgff$end) - pmax(segStart, sgff$start) + 1) / sgff$length
          
          whf = which(overlap>0)
          segScore[j, p("feature")] = paste(unique(sgff$Name[whf]), collapse=", ")
          
          ## fraction of overlap with features
          segScore[j, p("overlap")] = max(overlap, 0)
            
          ## Please see also man page:
          ## whs = all features that contain start of this segment:
          whs = which(segStart >= sgff$start & segStart < sgff$end)
          if(length(whs)>0) {
            segScore[j, p("dist.start2feat")] = min(segStart - sgff$start[whs])
          } else {
            segScore[j, p("dist.start2feat")] = posMin(segStart - sgff$end)
          }
          ## whe = all features that contain end of this segment:
          whe = which(segEnd >= sgff$start & segEnd < sgff$end)
          if(length(whe)>0) {
            segScore[j, p("dist.end2feat")] = min(sgff$end[whe] - segEnd)
          } else {
            segScore[j, p("dist.end2feat")] = posMin(sgff$start - segEnd)
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
