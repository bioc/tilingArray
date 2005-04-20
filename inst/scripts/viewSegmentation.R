options(error=recover, warn=2)

library("tilingArray")
library("arrayMagic")   ## for write.htmltable
library("geneplotter")  ## for savetiff

source("readSegments.R")
source("colorRamp.R") ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

rt = "tot"

outdir = file.path(indir[rt], "viz")

if(!file.exists(outdir) || !file.info(outdir)$isdir)
  stop(paste("Output directory", outdir, "does not exist."))

#### Generic plot
X11(); grid.newpage()
e = get(rt)
plotAlongChrom2(chr=3, coord = c(130, 150)*1e3, segRes = e, segScore = get("segScore", e), 
     gff = gff) # , nrBasesPerSeg = 2000)

stop()

## new transcripts:
sel = which( (segScore$same.feature=="") &
             (segScore$frac.dup < 0.2) &
     ##      !is.na(segScore$oppo.feature)
             (segScore$length > 200) & 
             (segScore$chr != 17) &
             (segScore$pt < 1e-10) )


## order by level
ord = order(segScore$level[sel], decreasing=TRUE)
sel = sel[ord]

outtab = cbind(plot=I(paste('<a href="', sel, '.tiff">', sel, '</a>', sep="")),
               segScore[sel, ])
  
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
  cat(i, ":  ", segScore$chr[s], segScore$strand[s], "  ",
      segScore$start[s], "...", segScore$end[s], "   ", sep="")

  tmpf = paste(tempfile(), "pdf", sep=".")
  pdf(file=tmpf, width=12, height=6)
  grid.newpage()
  plotAlongChrom2(chr=as.numeric(segScore$chr[s]),
                  coord=c(max(segScore$start[s]-2e4, 0),
                          segScore$end[s]+2e4),
                  highlight= list(coord=c(segScore$start[s], segScore$end[s]),
                                  strand=segScore$strand[s]),
                  segRes = segRes,
                  segScore = segScore, gff = gff)
  dev.off()
  cmd = paste("convert -density 120", tmpf, "-compress RLE",
              file.path(outdir, paste(s, "tiff", sep=".")), "; rm", tmpf, "&")
  system(cmd)
}
