options(error=recover)

library("tilingArray")
library("arrayMagic")   ## for write.htmltable
library("geneplotter")  ## for savetiff

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

outtab = cbind(plot=I(paste('<a href="', sel, '.tiff">', sel, '</a>', sep="")),
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

cat("Writing", length(sel), "alongChrom images.\n")
for(i in seq(along=sel)) {
  s = sel[i]
  cat(i, segScore$chr[s], segScore$strand[s], "  ",
      segScore$start[s], "...", segScore$end[s], "   ", sep="")

  tmpf = paste(tempfile(), "pdf", sep=".")
  pdf(file=tmpf, width=12, height=6)
  grid.newpage()
  plotAlongChrom2(chr=as.numeric(segScore$chr[s]),
                  coord=c(max(segScore$start[s]-2e4, 0),
                          segScore$end[s]+2e4),
                  highlight= list(coord=(segScore$start[s]+segScore$end[s])/2,
                                  strand=segScore$strand[s]),
                  segRes = segRes,
                  segScore = segScore, gff = gff)
  dev.off()
  cmd = paste("convert -density 120", tmpf, "-compress RLE",
              file.path(outdir, paste(s, "tiff", sep=".")), "; rm", tmpf, "&")
  system(cmd)
}
