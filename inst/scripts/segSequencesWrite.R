## write FASTA files with segments so they can be blasted against
## another genome
options(error=recover)

library("tilingArray")
library("matchprobes")

indir  = "segmentation-050209v4"
seqdir = "SGD"

outdir = file.path(indir, "fasta")

if(!file.exists(outdir) || !file.info(outdir)$isdir)
  stop(paste("Output directory", outdir, "does not exist."))

if(!exists("segScore"))
  load(file.path(indir, "segScore.rda"))

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

## We want to distnguish three groups: annotated transcripts,
## not annotated transcripts, and not annotated not transcribed
## sequences.

## At this point, we simply write out all segments:

cat("Writing", nrow(segScore), "sequences.\n")
con = file(file.path(outdir, "segments.fsa"), open="wt")


stopifnot(all(segScore$strand[s] %in% c("+", "-")))
for(s in 1:nrow(segScore)) {
  sequence = substr(fsa[[paste(segScore$chr[s])]],
                 start = segScore$start[s],
                 stop  = segScore$end[s])
  if(segScore$strand[s] == "-")
    sequence = complementSeq(sequence)
  
  cat(">", s, "\n", sequence, "\n", sep="", file=con, append=TRUE)
}
close(con)
