options(warn=0)
library("tilingArray")

nrBasesPerSeg = 1400

options(error=recover, warn=2)
if(!exists("probeAnno"))
  load("probeAnno.rda")

outdirList = c("polyA2" = "seg-polyA-050525",
  "tot"    = "seg-tot-050525",
  "tot2"   = "seg-tot2-050525")

chrstr = paste(rep(1:17, each=2),
               rep(c("+", "-"), 17), sep=".")

## ------------------------------------------------------------
## fuction subsample
## ------------------------------------------------------------
subSample = function(x, step=7) {
  stopifnot(all(diff(x)>=0))
  keep = rep(TRUE, length(x))
  done = FALSE

  i = 1
  while(!done) {
    xk =x[keep] 
    ## distance to second next neighbour
    d2 = xk[3:length(xk)] - xk[1:(length(xk)-2)]
    tooClose = c(FALSE, (d2<step), FALSE)
    ## but only if the preceding probe was not not too close
    ## itself
    tooClose = tooClose & !c(tooClose[-1], FALSE)
    if(any(tooClose)) {
      keep[keep][tooClose] = FALSE
    } else {
      done = TRUE
    }
    cat(i, sum(!keep), "\n")
  }
  browser()

  
  return(keep)
}

## ------------------------------------------------------------
## main
## ------------------------------------------------------------

for(outdir in outdirList) {
  if(!file.exists(outdir))
    stop(paste("Directory '", outdir, "' does not exist.", sep=""))
  if(!file.info(outdir)$isdir)
    stop(paste("'", outdir, "' must be a directory.", sep=""))

  for(ichr in seq(along=chrstr)) {
    chr = chrstr[ichr]
    #if(ichr==1) {
    #  load(file.path(outdir, "xn.rda"))
    #  lxj = exprs(xn)
    #  stopifnot(ncol(lxj) %in% c(2,3))
    #} 
    cat("OOOOOOPPPSSS!!!!\n")
    
    datfn = file.path(outdir, paste(chr, ".rda", sep=""))
    if(!file.exists(datfn)) {
      cat("working on", datfn, "\n")
      ## Lock it
      con = file(datfn, open="wt")
      writeLines(date(), con)
      close(con)
    
      ## get data
      sta = get(paste(chr, "start", sep="."), probeAnno)
      ord = order(sta)
      sta = sta[ord]
      end = get(paste(chr, "end", sep="."), probeAnno)[ord]
      ind = get(paste(chr, "index", sep="."), probeAnno)[ord]
      uni = get(paste(chr, "unique", sep="."), probeAnno)[ord]
      yraw = lxj[ind, ]

      ss = subSample(sta)
      
      ## use approx
      mids = (sta+end)/2
      dat  = list(start = sta, end = end, yraw = yraw, unique = uni,
        x=seq(min(mids), max(mids), by=8))
      dat$y        = apply(yraw, 2, function(y) { approx(mids, y, xout=dat$x, ties=mean)$y } )
      dat$xunique  = as.logical(round(approx(mids, uni, xout=dat$x, ties=mean)$y))
     
      ## see plotFeatSize: the longest structure CDS in yeast is
      ## 15k bases long, corresponding to about 2000 consecutive probes
      maxk  = 1500
      maxcp = round(diff(range(dat$x)) / nrBasesPerSeg)
      
      seg = findSegments(dat$y, maxk=maxk, maxcp=maxcp, verbose=99)
      
      save(seg, dat, file=datfn, compress=TRUE)
    } else {
      cat(datfn, "exists, skipping.\n")
    } ## if
  } # for chr
} ## for outdir

