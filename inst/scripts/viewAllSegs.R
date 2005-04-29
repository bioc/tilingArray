options(error=recover, warn=0)

library("tilingArray")

source("scripts/readSegments.R")
## source("colorRamp.R")  ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

rt = "polyA"
alongChromWidth = 25e3
alongChromStep  = 10e3
nrChr = 16

outdir = file.path(indir[rt], "viz")
if(!file.exists(outdir) || !file.info(outdir)$isdir)
  stop(paste("Output directory", outdir, "does not exist."))


## this function maps a chromosome number and start and end coordinates
## to a file name
mapCoord2Plot = function (chr, start, end) {
  mid  = (start+end-alongChromWidth)/2
  mid[mid<0]=0
  pst  = as.integer(alongChromStep/1e3*round(mid/alongChromStep))
  sprintf("%02d_%04d", as.integer(chr), pst)
}

##
## Write the GFF-table
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
write.table(outtab, file=file.path(outdir, "gff.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


##
## Write the pictures
## 

e = get(rt)
for(chr in 1:nrChr) {
  startPoints = seq(0, max(gff$end[gff$chr==chr]), by=alongChromStep)
  conv = character(length(startPoints))
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
    conv[i] = paste("convert -density 120", pdfname, "-quality 100", pixname)
  }

  scriptfn = sprintf("viewAllSeqs%02d.sh", as.integer(chr))
  conv = c("#!/bin/sh", conv)
  writeLines(conv, con=scriptfn)
  cmd  =  paste("chmod a+x ", scriptfn, "; ./", scriptfn, " &", sep="")
  cat("\n", cmd, "\n")
  system(cmd)
}
