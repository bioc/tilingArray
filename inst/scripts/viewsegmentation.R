options(error=recover, warn=2)

library("tilingArray")
source("colorRamp.R")
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

if(!exists("gff"))
  load("gff.rda")

indir = "segmentation-050209v4"

if(!exists("segRes")) {
  chrs = 1:17
  segRes  = new.env()
  cat("Loading ")
  for(chr in chrs) {
    for(strand in c("+", "-")) {
      fn = paste(chr, strand, "rda", sep=".")
      cat(fn, "")
      load(file.path(indir, fn))  
      assign(paste(chr, strand, "seg", sep="."), seg, envir=segRes)
      assign(paste(chr, strand, "dat", sep="."), dat, envir=segRes)
    }
  } ## for chr
  load(file.path(indir, "segScore.rda"))
} ## if
cat("\n")



#### Generic plot
##plotAlongChrom2(chr=1, coord = c(0, 230)*1e3, segRes = segRes,
##     gff = gff, nrBasesPerSeg = 1500)



## antisense segments:
sel = which(is.na(segScore$same.feature) & (segScore$frac.dup < .2)
         & !is.na(segScore$oppo.feature) & (segScore$chr != 17) )

## not annotated
##sel = which(is.na(segScore$same.feature))


ord = order(segScore$pt[sel])
sel = sel[ord]

for(s in sel) {
  cat(segScore$chr[s], segScore$strand[s], "  ",
      segScore$start[s], "...", segScore$end[s], "\n", sep="")
  grid.newpage()
  plotAlongChrom2(chr=as.numeric(segScore$chr[s]),
                  coord=c(max(segScore$start[s]-2e4, 0),
                          segScore$end[s]+2e4), 
                  segRes = segRes,
                  segScore = segScore, gff = gff)
  ## locator(n=1)
  readline("Type <enter>") 
}
