## write FASTA files with segments so they can be blasted against
## another genome
library("tilingArray")

interact = TRUE
rnaTypes  = c("seg-polyA-050525", "seg-tot-050525", "seg-tot2-050525")[1]
outfile = "tableSegments"

source(scriptsDir("readSegments.R"))
source(scriptsDir("categorizeSegments.R"))
source(scriptsDir("calcThreshold.R"))

options(error=recover)

if(!exists("gff"))
  load("probeAnno.rda")

segScoreFile = "segScore-1500.rda"
seqDir = "SGD"

outdir = "fasta"
what   = c("all","novel-isolated")[2]

for(rt in rnaTypes) {
  df = file.path(rt, outdir)
  if(!file.exists(df)) {
    cat("Creating output directory", df, ".\n")
    dir.create(df)
  }
  if(!file.info(df)$isdir)
    stop(paste(outdir, "is not a directory."))
}

if(!exists("fsa")) {
  fsa = new.env()
  fsa.files = paste("chr", c(sapply(1:16, function(n) sprintf("%02d", n)), "mt"),
    ".fsa", sep="")
  for(i in seq(along=fsa.files)) {
    s = readLines(file.path(seqDir, fsa.files[i]))
    s = paste(s[-1], collapse="")
    assign(paste(i), s, envir=fsa)
    cat(fsa.files[i], ": ", nchar(s), "\n", sep="")
  }
}


## Check if the sequence lengths found here coincide with the end of
## the telomere in the GFF table. If yes, all is well!
chrLengths = sapply(fsa, nchar)
chrLengths = chrLengths[order(as.numeric(names(chrLengths)))]

## double-check
sgff = gff[ gff$feature=="telomere", ]
for(i in 1:16) {
  w = which(sgff$chr==i)
  stopifnot(length(w)==2)
  stopifnot(chrLengths[i]==sgff$end[w[2]])
}

## write segment seqs:
for(rt in rnaTypes) {
  fout = file.path(rt, outdir, paste(what, "fsa", sep="."))

  segScore = categorizeSegments(get(rt))
  wh = switch(what,
    "all"            = {
      1:nrow(segScore)
    },
    "novel-isolated" = {
      which(segScore$category=="novel isolated - filtered")
    },
    stop("Zapperlot")
  )

  cat("Writing", length(wh), "sequences to", fout, ".\n")
  con = file(fout, open="wt")

  stopifnot(all(segScore$strand %in% c("+", "-")))
  for(s in wh) {
    sequence = substr(fsa[[paste(segScore$chr[s])]],
      start = segScore$start[s],
      stop  = segScore$end[s])
    cat(">", s, "\n", sequence, "\n", sep="", file=con, append=TRUE)
  }
  close(con)
}
