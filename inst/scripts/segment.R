options(warn=0)
library("prada")  ### for shorth
library("tilingArray")

nrBasesPerSeg = 1000

options(error=recover, warn=2)
if(!exists("probeAnno"))
  load("probeAnno.rda")

normalize = function(fn, x) {
  k    = which(x$File==fn);
  jref = which(x$Hybe %in% c(12,23,24))
  stopifnot(length(jref)==3, length(k)==1)
  return(log(exprs(x)[,k], 2) - rowMeans(log(exprs(x)[,jref], 2)))
}

if(!exists("lxj")) {
  load("x.Rdata")
  lxj = normalize("030505_totcDNA_15ug_affy.cel.gz", x)
  outdir = "segmentation-050305"
  ## lxj = normalize("050209_mRNAx4_30min_re-hybe_RH6.cel.gz", x)
  ## outdir = "segmentation-050209v4"
  rm(x)
  gc()
}

igvalues = lxj[probeAnno$probe$no_feature == "no"]
if(interactive())
  hist(igvalues, col="orange", breaks=100)
baseline = shorth(igvalues)

chrstr = paste(rep(1:17, each=2),
             rep(c("+", "-"), 17), sep=".")

## Output Directory
if(!file.exists(outdir))
  stop(paste("Directory '", outdir, "' does not exist.", sep=""))
if(!file.info(outdir)$isdir)
  stop(paste("'", outdir, "' must be a directory.", sep=""))

for(chr in chrstr) {
  datfn = file.path(outdir, paste(chr, ".rda", sep=""))
  if(!file.exists(datfn)) {
    cat("working on", datfn, "\n")
    ## Lock it
    con = file(datfn, open="wt")
    writeLines(date(), con)
    close(con)
    
    ## get data
    ind = get(paste(chr, "index", sep="."), probeAnno)
    sta = get(paste(chr, "start", sep="."), probeAnno)
    uni = get(paste(chr, "unique", sep="."), probeAnno)

    yraw = lxj[ind]-baseline
    
    ## average identical probes
    ## tm = tapply(yraw, sta, mean)
    ## mt = match(names(tm), paste(sta))
    ## stopifnot(!any(is.na(mt)))

    ## use approx
    dat     = approx(sta, yraw,   xout=seq(min(sta), max(sta), by=8), ties=mean)
    xunique = as.logical(round(approx(sta, uni, xout=dat$x, ties=mean)$y))

    dat = c(dat, list(xunique = xunique, start = sta, yraw = yraw, unique = uni,
                      baseline = baseline))
                 
    ## see plotFeatSize: the longest structure CDS in yeast is
    ## 15k bases long, corresponding to about 2000 consecutive probes
    maxk  = 1500
    maxcp = round(diff(range(dat$x)) / nrBasesPerSeg)

    seg = findsegments(dat$y, maxk=maxk, maxcp=maxcp, verbose=99)
      
    save(seg, dat, file=datfn)
  } else {
    cat(datfn, "exists, skipping.\n")
  } ## if
} # for chr


