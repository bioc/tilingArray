
## gff and chrSeqname into an environment or object?

plotAlongChrom = function(y, chr, from, to, extend=0, gff, chrSeqname,
  probeAnno, probeLength=25) {
  if(is.character(chr)) {
    chr = match(chr, chrSeqname)
    stopifnot(!is.na(chr))
  }
    
  if(from>to) {
    tmp=from;to=from;from=tmp
  }
  start  = from-extend
  end    = to+extend

  rgy  = quantile(y, probs=c(0.01, 0.99))
  y    = y - rgy[1]
  maxy = diff(rgy)

  thegff = gff[ (gff$seqname == chrSeqname[chr]) &
                (gff$feature == "CDS"), ]
  stopifnot(!is.null(gff$Name))
  
  plot(c(start, end), c(-1,1)*maxy, type="n", 
       xlab=paste("Chr", chr, sep=""), ylab="")
  abline(h=0, col="grey")

  ## Plot expression data
  ## Watson strand = antisense = "+" = typically plotted above
  ## Crick strand = sense = "-" = typically plotted below
  
  strandnames = c("-", "+")
  colors = c("-" = "#1F78B4", "+" = "#33A02C")
  sgn    = c("-" = -1, "+" = +1)
  stopifnot(identical(strandnames, names(colors)),
            identical(strandnames, names(sgn)),
            all(as.character(thegff$strand) %in% names(sgn)))
  
  for(strand in strandnames) {
    pos = get(paste(chr, strand, "start", sep="."), envir=probeAnno)
    ind = get(paste(chr, strand, "index", sep="."), envir=probeAnno)
    uni = get(paste(chr, strand, "unique", sep="."), envir=probeAnno)
    
    probesel = which(((pos >= start) & (pos+probeLength <= end)) |
                ((pos+probeLength <= start) & (pos >= end)))

    points(pos[probesel], y[ind[probesel]] * sgn[strand], pch=20, cex=0.5,
           col= ifelse(uni[probesel], colors[strand], "grey"))
  }

  ## Plot features
  featsel = logical(nrow(thegff))
  for(i in 1:nrow(thegff)) {
    ig = seq(thegff$start[i], thegff$end[i])
    featsel[i] = any(ig >= start & ig <= end)
  }
  if(any(featsel)) {
    strd = as.character(thegff$strand[featsel])
    y  = sgn[strd] * 0.2 * maxy
    ax = cbind(thegff$start[featsel], thegff$end[featsel])
    ax[strd=="-", ] = ax[strd=="-", 2:1]
    arrows(x0=ax[,1], y0=y, x1=ax[,2], y1=y, lwd=2, length=0.1, col="black", code=2)
        ##   col=colors[strd]
    text(x=(thegff$start[featsel]+thegff$end[featsel])/2, y=-1+seq(along=which(featsel))%%3,
         labels=thegff$Name[featsel])
  }
}

