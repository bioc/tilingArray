library("tilingArray")

graphics.off()
options(error=recover, warn=2)

rnaTypes = c("seg-polyA-050811", "seg-tot-050811")
outfile  = "lowTranscribedRegions"

source(scriptsDir("readSegments.R")) 

if(!interact){
  sink(paste(outfile, ".txt", sep=""))
  cat("Made on", date(), "\n\n")
}


if(!interact)
  sink()
