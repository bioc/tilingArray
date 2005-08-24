options(warn=0)
library("tilingArray")

nrBasesPerSeg = 1400

options(error=recover, warn=2)
if(!exists("probeAnno"))
  load("probeAnno.rda")

##outdirs = c("seg-polyA-050525", "seg-tot-050525", "seg-dir-050721" ,
##  "seg-odT-050801", "seg-polyA-050804")[3:5]

## these use the same normalized data as the ones above, but they
## use newer probeAnno file (0508)
outdirs = c("seg-polyA-050811", "seg-tot-050811", "seg-dir-050811",
  "seg-odT-050811", "seg-polyA0420-050811")

## these use the same normalized data as the ones above, but they
## use newer probeAnno file (0508) and the new findSegments (tilingArray
outdirs = c("seg-polyA-050824", "seg-tot-050824", "seg-dir-050824",
  "seg-odT-050824", "seg-polyA0420-050824")

##  for(d in outdirs) {
##    dir.create(d); dir.create(file.path(d, "viz")) }

chrstr = paste(rep(1:17, each=2),
               rep(c("+", "-"), 17), sep=".")

allPM = unique(unlist(lapply(chrstr, function(chr)
  get(paste(chr, "index", sep="."), probeAnno))))

## ------------------------------------------------------------
## main
## ------------------------------------------------------------

for(outdir in outdirs) {
  if(!file.exists(outdir))
    stop(paste("Directory '", outdir, "' does not exist.", sep=""))
  if(!file.info(outdir)$isdir)
    stop(paste("'", outdir, "' must be a directory.", sep=""))

  for(ichr in seq(along=chrstr)) {
    chr = chrstr[ichr]
    if(ichr==1) {
      load(file.path(outdir, "xn.rda"))
      lxj = exprs(xn)
      refSigThresh = quantile(refSig[allPM], probs=0.05)
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
      end = get(paste(chr, "end", sep="."), probeAnno)
      mid = (sta+end)/2
      ord = order(mid)

      sta = sta[ord]
      end = end[ord]
      mid = mid[ord]
      ind = get(paste(chr, "index", sep="."), probeAnno)[ord]
      uni = get(paste(chr, "unique", sep="."), probeAnno)[ord]
      y   = lxj[ind,, drop=FALSE]

      ## eliminate the bad probes
      sel = (refSig[ind] >= refSigThresh)

      ## subsample to keep spacing around 8
      ss = sampleStep( mid[sel], step=7)
      
      dat  = list(start  = sta[sel],
                  end    = end[sel],
                  y      = y[sel,, drop=FALSE],
                  unique = uni[sel],
                  ss     = ss)
     
      ## see plotFeatSize: the longest structure CDS in yeast is
      ## 15k bases long, corresponding to about 2000 consecutive probes
      maxk  = 1500
      maxcp = round( (max(end)-min(sta)) / nrBasesPerSeg)
      
      seg = findSegments(dat$y[ss,,drop=FALSE], maxk=maxk, maxcp=maxcp, verbose=99)
      
      save(seg, dat, file=datfn, compress=TRUE)
    } else {
      cat(datfn, "exists, skipping.\n")
    } ## if
  } # for chr
} ## for outdir

