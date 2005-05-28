##
## Note: to produce the input files (.fasta) for this script, run
## segSequencesWrite.R
##
## S. paradoxus, S. mikatae, S. bayanus

segdir = file.path(c("seg-polyA-050525", "seg-tot-050525", "seg-tot2-050525"), "fasta")


theDir   = file.path("OtherSpecies", c("S.bayanus", "S.mikatae", "S.paradoxus"))
theFiles = c("Sbay_contigs.fasta", "Smik_contigs.fasta", "Spar_contigs.fasta")


cat("Please run:\n")
for(i in seq(along=theDir)) {
  cat("cd", theDir[i], "\n",
      "formatdb -p F -i", theFiles[i], "-o T\n",
      "cd -\n\n")
}


bsubfile = "blast_threefunghi.sh"
bsubcon  = file(bsubfile, open="wt")

for(j in seq(along=segdir)){
  for(i in seq(along=theDir)) {
    jobfile = sprintf("blast_threefunghi_%d_%d.sh", as.integer(j), as.integer(i))
    con = file(jobfile, open="wt")
    cat("blastall -p blastn -d ", file.path(theDir[i], theFiles[i]),
      " -i ", file.path(segdir[j], "segments.fsa"), " -m 8 -e 1e-10 -a 4 -o ", 
      file.path(segdir[j], sub(".fasta", ".out", theFiles[i])), "\n\n", sep="", file=con)
    close(con)
    system(paste("chmod 755", jobfile))
    cat("bsub -q research ./", jobfile, "\n", sep="", file=bsubcon)
  }
}

close(bsubcon)
system(paste("chmod 755", bsubfile))
