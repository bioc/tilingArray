options(error=recover, warn=2)

library("tilingArray")
library("arrayMagic")

source("colorRamp.R") ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")
source("~/madman/Rpacks/arrayMagic/R/write.htmltable.R")

if(!exists("gff"))
  load("gff.rda")

indir  = "segmentation-050209v4"
outdir = file.path(indir, "viz")

if(!file.exists(outdir) || !file.info(outdir)$isdir)
  stop(paste("Output directory", outdir, "does not exist."))

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
sel = which(is.na(segScore$same.feature) & (segScore$frac.dup < 0.2) 
     ##    & !is.na(segScore$oppo.feature)
           & (segScore$end-segScore$start > 200) 
           & (segScore$chr != 17) & (segScore$pt < 1e-10))


## order by level
ord = order(segScore$level[sel], decreasing=TRUE)
sel = sel[ord]

sel = sel[1:20]

outtab = cbind(plot=I(paste('<a href="', sel, '.png">plot</a>', sep="")),
               segScore[sel, 1:4],
               length=as.integer(segScore$end[sel]-segScore$start[sel]+1),
               segScore[sel, 5:ncol(segScore)])
  
for(i in which(colnames(outtab) %in% c("start", "end",
     "same.dist.start2feat", "same.dist.end2feat",
     "oppo.dist.start2feat", "oppo.dist.end2feat")))
  outtab[[i]] = as.integer(outtab[[i]])

outtab = outtab[-which(colnames(outtab) %in% c("same.feature", "same.overlap"))]

write.htmltable(outtab,
   filename=file.path(outdir, "newtranscripts"),
   title=paste(length(sel), " segments with candidates for new transcripts (", indir, ")", sep=""))

for(s in sel) {
  cat(segScore$chr[s], segScore$strand[s], "  ",
      segScore$start[s], "...", segScore$end[s], "\n", sep="")
  png(file=file.path(outdir, paste(s, "png", sep=".")), width=1280, height=640)
  grid.newpage()
  plotAlongChrom2(chr=as.numeric(segScore$chr[s]),
                  coord=c(max(segScore$start[s]-2e4, 0),
                          segScore$end[s]+2e4),
                  highlight= list(coord=(segScore$start[s]+segScore$end[s])/2,
                                  strand=segScore$strand[s]),
                  segRes = segRes,
                  segScore = segScore, gff = gff)
  dev.off()
  
  ## locator(n=1)
  ## readline("Type <enter>") 
}
