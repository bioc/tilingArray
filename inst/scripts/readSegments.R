if(!exists("gff")) {
  cat("Loading probeAnno.rda\n")
  load("probeAnno.rda")
  gff$Name = getAttributeField(gff$attributes, "Name")
  theID    = getAttributeField(gff$attributes, "ID")
  stopifnot(all(gff$Name == theID, na.rm=TRUE))
  gff$orf_classification = getAttributeField(gff$attributes, "orf_classification")
  gff$gene               = getAttributeField(gff$attributes, "gene")
}

rnaTypes = c("polyA", "polyA2", "tot")[2:3]
longNames=c(polyA="poly-A RNA single enriched", polyA2="poly-A RNA", tot="total RNA")

indir = c("polyA"="segmentation-3polyA", "polyA2"="seg-polyA-050428",
  "tot"="seg-tot-050421")

for(rt in rnaTypes) {
  if(!exists(rt)) {
    assign(rt, new.env())
    chrs = 1:17
    cat("Loading", rt, ": ")
    for(chr in chrs) {
      cat(chr, "")
      for(strand in c("+", "-")) {
        fn = paste(chr, strand, "rda", sep=".")
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


