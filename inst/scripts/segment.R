options(warn=0)
library("prada")  ### for shorth
library("tilingArray")

nrBasesPerSeg = 1000

options(error=recover, warn=2)
if(!exists("probeAnno"))
  load("probeAnno.rda")

if(!exists("x")) 
  load("x.rda")

what = c("polyA", "tot")[1]

switch(what,
  "polyA" = {
    fn = c("041203_S96_polyAx1_RH6.cel.gz",
      "050209_mRNAx4_30min_re-hybe_RH6.cel.gz",
      "050218_polyA-RNA_RH6_4x15min.cel.gz")
    outdir = "seg-polyA-050418"
  },
  "tot"    = {
    fn = c("050409_totcDNA_14ug_no52.cel.gz",
      "030505_totcDNA_15ug_affy.cel.gz",
      "050415_totcDNA_20ug_Affy11.cel.gz")
    outdir = "seg-tot-050418"
  },
  stop(paste("Bummer:", what))
)

lxj = exprs(x)[, fn]
stopifnot(ncol(lxj)==3)

hybeType = "Reverse"

igvalues = lxj[get(paste("probe", hybeType, sep=""), envir=probeAnno)$no_feature == "no", ]
baseline = shorth(igvalues, na.rm=TRUE)
if(interactive())
  hist(igvalues, col="orange", breaks=100, main=sprintf("baseline=%3.1f", baseline))

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

    yraw = lxj[ind, ]-baseline

    ## use approx
    dat = list(x=seq(min(sta), max(sta), by=8))
    dat$y        = apply(yraw, 2, function(y) { approx(sta, y, xout=dat$x, ties=mean)$y } )
    dat$xunique  = as.logical(round(approx(sta, uni, xout=dat$x, ties=mean)$y))
    dat$start    = sta
    dat$yraw     = yraw
    dat$unique   = uni
    dat$baseline = baseline
     
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


