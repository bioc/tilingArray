## Look at the fraction of duplicated probes of segments


library("tilingArray")

indir = "segmentation-050209v4"
if(!exists("segScore"))
  load(file.path(indir, "segScore.rda"))

hist(segScore$frac.dup, col="orange", breaks=seq(0,1,by=0.02), 
     main=paste("fraction of duplicated probes in segments\n", indir))

dev.copy(pdf, file="segmentDuplication.pdf")
dev.off()
