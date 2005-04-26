## S. paradoxus, S. mikatae, S. bayanus

segdir = "segmentation-3polyA/fasta"


theDir   = file.path("OtherSpecies", c("S.bayanus", "S.mikatae", "S.paradoxus"))
theFiles = c("Sbay_contigs.fasta", "Smik_contigs.fasta", "Spar_contigs.fasta")

for(i in seq(along=theDir)) {
  cat("Please run:\n",
      "cd", theDir[i], "\n",
      "formatdb -p F -i", theFiles[i], "-o T\n",
      "cd -\n\n")
}
   
for(i in seq(along=theDir)) {
  cat("blastall -p blastn -d ", file.path(theDir[i], theFiles[i]),
      " -i ", file.path(segdir, "segments.fsa"), " -m 8 -e 1e-10 -a 4 -o ", 
      file.path(segdir, sub(".fasta", ".out", theFiles[i])), "\n\n", sep="")
}