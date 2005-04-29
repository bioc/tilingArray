options(error=recover, warn=0)

library("tilingArray")

source("scripts/readSegments.R")
## source("colorRamp.R")  ## can go with R 2.1
## source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

rt = "polyA"

#### Generic plot
graphics.off(); X11(width=15, height=8); grid.newpage()
e = get(rt)

stopifnot("Name" %in% names(gff))
w = which(gff$Name=="YAL003W" & gff$feature=="gene")
stopifnot(length(w)==1)

if(TRUE) {
  ## with segRes environment
  plotAlongChrom2(which(gff$seqname[w]==chrSeqname), coord = c(gff$start[w]-1e4, gff$end[w]+1e4),
                nrBasesPerSeg=1500, segRes = e,
                ## segScore = get("segScore", e), 
                gff = gff, highlight= list(coord=c(142621, 143365),strand="+"))
} else {
  if(!exists("a"))load("a.rda")
  if(!exists("probeAnno"))load("probeAnno.rda")
  plotAlongChrom2(which(gff$seqname[w]==chrSeqname), coord = c(gff$start[w]-1e4, gff$end[w]+1e4),
                  nrBasesPerSeg=1000, y = log(exprs(a)[, "05242_totRNA_15ugS96_dir#3.cel.gz"], 2),
                  probeAnno = probeAnno,
                  gff = gff, highlight= list(coord=c(142621, 143365),strand="+"))
}
