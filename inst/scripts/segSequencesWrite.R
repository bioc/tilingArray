## write FASTA files with segments so they can be blasted against
## another genome

options(error=recover)

library("tilingArray")
## library("matchprobes")

if(!exists("gff"))
  load("probeAnno.rda")

indir  = "segmentation-3polyA"
segScoreFile = "segScore-1500.rda"
seqdir = "SGD"

outdir = file.path(indir, "fasta")

if(!file.exists(outdir) || !file.info(outdir)$isdir)
  stop(paste("Output directory", outdir, "does not exist."))

if(!exists("segScore"))
  load(file.path(indir, segScoreFile))

if(!exists("fsa")) {
  fsa = new.env()
  fsa.files = paste("chr", c(sapply(1:16, function(n) sprintf("%02d", n)), "mt"),
    ".fsa", sep="")
  for(i in seq(along=fsa.files)) {
    s = readLines(file.path(seqdir, fsa.files[i]))
    s = paste(s[-1], collapse="")
    assign(paste(i), s, envir=fsa)
    cat(fsa.files[i], ": ", nchar(s), "\n", sep="")
  }
}


### Check if the sequence lengths found here coincide with the end of
## the telomere in the GFF table. If yes, all is well!

chrLengths = sapply(fsa, nchar)
chrLengths = chrLengths[order(as.numeric(names(chrLengths)))]

## double-check
sgff = gff[ gff$feature=="telomere", ]
for(i in 1:16) {
  w = which(sgff$seqname==chrSeqname[i])
  stopifnot(length(w)==2)
  stopifnot(chrLengths[i]==sgff$end[w[2]])
}

## We want to distinguish three groups: annotated transcripts,
## not annotated transcripts, and not annotated not transcribed
## sequences.
## At this point, we simply write out all segments:

cat("Writing", nrow(segScore), "sequences.\n")
con = file(file.path(outdir, "segments.fsa"), open="wt")


stopifnot(all(segScore$strand %in% c("+", "-")))
for(s in 1:nrow(segScore)) {
  sequence = substr(fsa[[paste(segScore$chr[s])]],
                 start = segScore$start[s],
                 stop  = segScore$end[s])
  cat(">", s, "\n", sequence, "\n", sep="", file=con, append=TRUE)
}
close(con)
