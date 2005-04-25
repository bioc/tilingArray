source("Figures/calcThreshold.R") 

if(!exists("gff")) {
  load("probeAnno.rda")
  gff$Name = getAttributeField(gff$attributes, "Name")
  theID   = getAttributeField(gff$attributes, "ID")
  stopifnot(all(gff$Name == theID, na.rm=TRUE))
  gff$orf_classification = getAttributeField(gff$attributes, "orf_classification")
  gff$gene = getAttributeField(gff$attributes, "gene")
}

rnaTypes = c("polyA", "tot")
indir = c("segmentation-3polyA", "seg-tot-050421")
names(indir) = rnaTypes

for(rt in rnaTypes) {
  if(!exists(rt)) {
    assign(rt, new.env())
    chrs = 1:16
    cat("Loading ", rt, ": ")
    for(chr in chrs) {
      for(strand in c("+", "-")) {
        fn = paste(chr, strand, "rda", sep=".")
        cat(fn, "")
        load(file.path(indir[rt], fn))  
        assign(paste(chr, strand, "seg", sep="."), seg, envir=get(rt))
        assign(paste(chr, strand, "dat", sep="."), dat, envir=get(rt))
      }
    } ## for chr
    fn="segScore-1500.rda"
    cat(fn, "\n")
    load(file.path(indir[rt], fn), envir=get(rt))
  } ## if
} ## rnaTypes

