options(error=recover, warn=0)

library("tilingArray")
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

## source("colorRamp.R")  ## can go with R 2.1
source("scripts/readSegments.R")
source("scripts/writeSegmentTable.R")

nrChr = 16

##
## Prepare the GFF-table
##
if(!"chr" %in% names(gff))
  gff$chr = match(gff$seqname, chrSeqname)
if(!"Ontology_term" %in% names(gff))
  gff$Ontology_term = getAttributeField(gff$attributes, "Ontology_term")
if(!"Note" %in% names(gff)) {
  tmp = getAttributeField(gff$attributes, "Note")
  tmp = gsub("%20", " ", tmp)
  tmp = gsub("%2C", ",", tmp)
  tmp = gsub("%3B", ";", tmp)
  gff$Note = tmp
}

outtab = gff[ gff$feature=="gene", c("chr", "start", "end",
                "strand", "Name", "gene", "orf_classification",
                "Note") ]
outtab$plotfile = mapCoord2Plot(outtab$chr, outtab$start, outtab$end)



##
## Loop over the RNA Types
##
convCmd = "#!/bin/sh"

for(rt in rnaTypes) {
  cat(">>>", rt, "<<<\n")
  
  ##
  ## Write the GFF-table
  ##
  outdir = file.path(indir[rt], "viz")
  if(!file.exists(outdir) || !file.info(outdir)$isdir)
    stop(paste("Output directory", outdir, "does not exist."))

  write.table(outtab, file=file.path(outdir, "gff.txt"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  ##
  ## Write the pictures
  ## 
  e = get(rt)
  for(chr in 1:nrChr) {
    startPoints = seq(0, max(gff$end[gff$chr==chr]), by=alongChromStep)
    for(i in seq(along=startPoints)) {
      start = startPoints[i]
      nm = mapCoord2Plot(chr, start+alongChromWidth/2, start+alongChromWidth/2)
      cat(chr, ":", as.integer(start), " ", sep="")
      pdfname = file.path(outdir, paste(nm, ".pdf", sep=""))
      pixname = file.path(outdir, paste(nm, ".jpg", sep=""))
      
      pdf(file=pdfname, width=12, height=6)
      grid.newpage()
      plotAlongChrom2(chr=chr, coord=c(start, start+alongChromWidth),
                      segRes   = e,
                      segScore = get("segScore", e), 
                      gff      = gff)
      dev.off()	
      convCmd = c(convCmd, paste("convert -density 120", pdfname, "-quality 100", pixname))
    }
  }
}


    
writeLines(convCmd, con="viewAllSeqs.sh")
system("chmod a+x viewAllSeqs.sh")
cat("\n\nYou can now run", viewAllSeqs.sh, "\n\n")

