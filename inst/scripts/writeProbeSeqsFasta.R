## 9.1.2005
## Write probe sequences into a FASTA format. Sequences are identified by their index
options(warn=2)
options(error=recover)

library("Scerevisiaetilingprobe")
if(!exists("idx")) {
  idx = as.integer(Scerevisiaetilingprobe$y*2560 + Scerevisiaetilingprobe$x + 1)
  stopifnot(identical(idx, 1:length(Scerevisiaetilingprobe$x)))
}

writeFasta = function(s, fname) {
  sid = sapply(idx[s], function(i) sprintf("%07d", i))
  seq = reverseSeq(Scerevisiaetilingprobe$sequence[s])
  pn = paste(">", sid, "\n", seq, sep="")
  writeLines(pn, con=fname)
}

sel  = which(!is.na(Scerevisiaetilingprobe$sequence))

ssel = sel[Scerevisiaetilingprobe$qualifier[sel] %in% c("CHR6_at", "CHR6_st") &
           Scerevisiaetilingprobe$destype[sel]   %in% c("-111", "111") ]

ssel = ssel[order(Scerevisiaetilingprobe$expos[ssel])]
ssel = ssel[(1:10)+length(ssel)/2]

jcut = round(seq(0, length(sel), length=11))
for(i in 2:length(jcut)){
  fn = paste("Scerevisiaetilingprobe", i, "fsa", sep=".")
  j1 = jcut[i-1]+1
  j2 = jcut[i]
  cat(fn, " ", j1, ":", j2, "\n", sep="")
  writeFasta(s = sel[j1:j2], fname=fn)
  cmd = paste("#!/bin/sh\n/usr/local/package/bin/blastall -p blastn -d SCgenome/SCgenome -i", fn,
            "-m 8 -e 0.01 -a 4 -o", gsub(".fsa", ".out", fn))
  writeLines(cmd, con = paste("Xtetris3d-V", i, sep=""))

  ## needs "force=TRUE" and "current working directory"
}

## writeFasta(s= ssel, "test.fsa")
## writeFasta(s=  sel, "Scerevisiaetilingprobe.fsa")
