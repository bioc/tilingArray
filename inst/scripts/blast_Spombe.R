
theDir  = "S.pombe"
infiles = paste("chromosome", 1:3, ".contig.embl", sep="")
outfile = "all.fsa"

con = file(file.path(theDir, outfile), open="wt")
for(i in seq(along=infiles)) {
  s      = readLines(file.path(theDir, infiles[i]))
  is.seq = (substr(s, start=1, stop=4)=="    ")
  ss     = s[is.seq]
  ss     = paste(substr(ss,  6, 15), substr(ss, 17, 26), substr(ss, 28, 37),
                 substr(ss, 39, 48), substr(ss, 50, 59), substr(ss, 61, 70), sep="") 
  writeLines(ss, con=con)
}
close(con)

segdir = "segmentation-050209v4/fasta/"

cat("Now you can run:\n")
cat("cd S.pombe\nformatdb -p F -i all.fsa -o T\n\n")
cat("cd ..\n")
cat("blastall -p blastn -d S.pombe/all.fsa -i ", segdir, "segments.fsa -m 8 -e 1e-10 -a 4 -o ", segdir, "segments.out\n", sep="")
