
## S. paradoxus, S. mikatae, S. bayanus

segdir = c("segmentation-3polyA/fasta",
  "seg-tot-050421/fasta")[2]


theDir   = file.path("OtherSpecies", c("S.bayanus", "S.mikatae", "S.paradoxus", "S.pombe"))
theFiles = c("Sbay_contigs.fasta", "Smik_contigs.fasta", "Spar_contigs.fasta", "Spom_all.fasta")


cat("Please run:\n")
for(i in seq(along=theDir)) {
  cat("cd", theDir[i], "\n",
      "formatdb -p F -i", theFiles[i], "-o T\n",
      "cd -\n\n")
}
   
for(i in seq(along=theDir)) {
  cat("blastall -p blastn -d ", file.path(theDir[i], theFiles[i]),
      " -i ", file.path(segdir, "segments.fsa"), " -m 8 -e 1e-10 -a 4 -o ", 
      file.path(segdir, sub(".fasta", ".out", theFiles[i])), "\n\n", sep="")
}
