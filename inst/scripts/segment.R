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
## main
## ------------------------------------------------------------

for(outdir in outdirList) {
  if(!file.exists(outdir))
    stop(paste("Directory '", outdir, "' does not exist.", sep=""))
  if(!file.info(outdir)$isdir)
    stop(paste("'", outdir, "' must be a directory.", sep=""))

  for(ichr in seq(along=chrstr)) {
    chr = chrstr[ichr]
    if(ichr==1) {
      load(file.path(outdir, "xn.rda"))
      lxj = exprs(xn)
      stopifnot(ncol(lxj) %in% c(2,3))
    } 
    
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
      y   = lxj[ind, ]

      ss = sampleStep( (sta+end)/2, step=7)
      
      dat  = list(start = sta, end = end, y = y, unique = uni, ss = ss)
     
      ## see plotFeatSize: the longest structure CDS in yeast is
      ## 15k bases long, corresponding to about 2000 consecutive probes
      maxk  = 1500
      maxcp = round( (max(end)-min(sta)) / nrBasesPerSeg)
      
      seg = findSegments(dat$y[ss, ], maxk=maxk, maxcp=maxcp, verbose=99)
      
      save(seg, dat, file=datfn, compress=TRUE)
    } else {
      cat(datfn, "exists, skipping.\n")
    } ## if
  } # for chr
} ## for outdir

