options(error=recover)

library("tilingArray")
library("matchprobes")
library("Biostrings")

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

  ## Check if the sequence lengths found here coincide with the end of
  ## the telomere in the GFF table. If yes, all is well!
  if(!exists("gff")) load("gff.rda")

  chrLengths = sapply(fsa, nchar)
  chrLengths = chrLengths[order(as.numeric(names(chrLengths)))]
  
  sgff = gff[ gff$feature=="telomere", ]
  for(i in 1:16) {
    w = which(sgff$seqname==chrSeqname[i])
    stopifnot(length(w)==2)
    cat(i, chrLengths[i] -  sgff$end[w[2]], "\n")
  }
}

segScore$orfcoverage = numeric(nrow(segScore))
for(chr in ls(fsa)) {
  sequence = get(paste(chr), fsa)
  for(strand in c("+", "-")) {
    a = NucleotideString(sequence, type="DNA", alphabet = DNAPatternAlphabet())
    if(strand=="-")
      a = reverseComplement(a)
      
    stopcodon = c(as.matrix(matchDNAPattern("TGA", a))[,1], 
      as.matrix(matchDNAPattern("TAG", a))[,1], 
      as.matrix(matchDNAPattern("TAA", a))[,1])
    startcodon = as.matrix(matchDNAPattern("ATG", a))[,1] 
    ## for each startcodon, find the next in frame stopcodom
    potentialORF = logical(nchar(a))
    for(s in startcodon) {
      d = (stopcodon-s)
      d = min(d[(d%%3==0) & d>0])
      if(is.finite(d))
        potentialORF[s + (0:(d+2))] = TRUE

    }

    ss    = segScore[ segScore$chr==chr & segScore$strand==strand, ]
    cover = mapply(function(i1, i2)mean(potentialORF[i1:i2]),
                   ss$start, ss$end)
    
    plot(ss$level, cover, col=c("black", "red")[1 + (ss$same.feature!="")],
         main=paste(chr, strand))
  }
}




