options(error=recover, warn=0)
interact =  !TRUE
writeGFF =  TRUE
library("tilingArray")

source("setScriptsDir.R")
source(functionsDir("plotAlongChrom.R"))


## rnaTypes  = c("seg-polyA-050811", "seg-tot-050811",
##  "seg-dir-050811" , "seg-odT-050811", "seg-polyA0420-050811")

rnaTypes  = c("seg-polyA-050909", "seg-tot-050909",
  "seg-dir-050909" , "seg-odT-050909", "seg-polyA0420-050909")

isDirect  = c(FALSE, FALSE, TRUE, FALSE, FALSE)

names(isDirect) = rnaTypes
setwd("/ebi/research/huber/Projects/tilingArray")

# check if rnaTypes were specified correctly:
stopifnot(length(isDirect)==length(rnaTypes),
          all(file.info(rnaTypes)$isdir),
          all(file.info(file.path(rnaTypes,"viz"))$isdir))  

source(scriptsDir("readSegments.R"))
source(scriptsDir("calcThreshold.R"))
source(scriptsDir("writeSegmentTable.R"))


nrChr = 16

##
## Prepare and write the GFF-table
##

if(!"Ontology_term" %in% names(gff))
  gff$Ontology_term = getAttributeField(gff$attributes, "Ontology_term")

if(!"Note" %in% names(gff)) {
  tmp = getAttributeField(gff$attributes, "Note")
  tmp = gsub("%20", " ", tmp)
  tmp = gsub("%2C", ",", tmp)
  tmp = gsub("%3B", ";", tmp)
  gff$Note = tmp
}

if (writeGFF) {
  outtab = gff[ gff$feature=="gene", c("chr", "start", "end", "strand",
                  "Name", "gene", "orf_classification","Note") ]
  outtab$plotfile = mapCoord2Plot(outtab$chr, outtab$start, outtab$end)
  write.table(outtab, file="gff.txt",
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

##
## Write the pictures
## 
for(rt in rnaTypes) {
  cat(">>>", rt, "<<<\n")
  convCmd = "#!/bin/sh"

  outdir = file.path(indir[rt], "viz")
  if(!file.exists(outdir) || !file.info(outdir)$isdir)
    stop(paste("Output directory", outdir, "does not exist."))
    
  allY = unlist(lapply(1:nrChr, function(chr) {
    lapply(c("-", "+"), function(strand) {
      dat = get(paste(chr, strand, "dat", sep="."), get(rt))
      dat[["y"]][ dat[["ss"]], ]
    })
  }))
  ylim = quantile(allY, probs=c(0.001, 0.999))
  
  for(chr in 1:nrChr) {
    startPoints = seq(0, max(gff[gff[, "chr"]==chr, "end"]), by=alongChromStep)
    for(i in seq(along=startPoints)) {
      start = startPoints[i]
      nm = mapCoord2Plot(chr, start+alongChromWidth/2, start+alongChromWidth/2)
      cat(chr, ":", as.integer(start), " ", sep="")
      pdfname = file.path(outdir, paste(nm, ".pdf", sep=""))
      pixname = file.path(outdir, paste(nm, ".jpg", sep=""))

      if(interact){
        x11(width=10, height=5.5)
      } else {
        pdf(file=pdfname, width=10, height=5.5)
      }
      grid.newpage()
      plotAlongChrom(chr=chr, coord=c(start, start+alongChromWidth),
                     segObj=get(rt), gff = gff, ylim=ylim, isDirect=isDirect[rt])
      if(!interact)
        dev.off()
      
      convCmd = c(convCmd, paste("convert -density 120", pdfname, "-quality 100", pixname))
    }
  }
  convCmd = c(convCmd, paste("chmod 664 ", outdir, "/*", sep=""))
  shellFile = paste("v", rt, "sh", sep=".")
  writeLines(convCmd, con=shellFile)
  system(paste("chmod a+x", shellFile))
  cat("\n\nYou can now run", shellFile, "\n\n")
}


