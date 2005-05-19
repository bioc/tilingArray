options(error=recover, warn=0)

library("tilingArray")

## source("colorRamp.R")  ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")


#### Generic plot
## graphics.off();
X11(width=15, height=8); grid.newpage()

##stopifnot("Name" %in% names(gff))
##w = which(gff$Name=="YAL003W" & gff$feature=="gene")
##stopifnot(length(w)==1)

if(!TRUE) {
  ## with segRes environment
  source("scripts/readSegments.R")
  rt = "polyA"
  plotAlongChrom2(which(gff$seqname[w]==chrSeqname), coord = c(gff$start[w]-1e4, gff$end[w]+1e4),
                nrBasesPerSeg=1500, segRes = get(rt),
                ## segScore = get("segScore", e), 
                gff = gff, highlight= list(coord=c(142621, 143365),strand="+"))
} else {
  ## if(!exists("a"))load("a.rda")
  if(!exists("xn"))load("seg-polyA-050518/xn.rda")
  if(!exists("probeAnno"))load("probeAnno.rda")
  fn = c("05_04_27_2xpolyA_NAP3.cel.gz",
    "05_04_26_2xpolyA_NAP2.cel.gz",
    "05_04_20_2xpolyA_NAP_2to1.cel.gz")
  ## y = log(exprs(a)[, fn[1]], 2)
  ## y = log(exprs(a)[, 1:3], 2)
  y = exprs(xn)
  plotAlongChrom2(## 6, coord = c(19000, 19200),
                  ## 7, coord = c(971000, 980000),
                  1, coord = c(184, 185)*1e3,
                  ## 3, coord = c(35, 50)*1e3,
                  y = y, 
                  probeAnno = probeAnno,
                  gff = gff)#  highlight= list(coord=c(142621, 143365),strand="+"))
}

##1160    6      +  19073  19089     16 4.211916                                 
##6682    7      + 973989 974021     32 4.714935                                 
##18512   9      + 290333 290365     32 4.556339                                 
##41610  11      + 613389 613445     56 3.857509                                 
