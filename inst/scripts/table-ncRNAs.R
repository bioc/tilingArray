if(!exists("gff")) {
  cat("Loading probeAnno.rda\n")
  load("probeAnno.rda")
}

f = c("ncRNA", "nc_primary_transcript",
      "rRNA", "snRNA",
      "snoRNA", "tRNA")

nm = unique(gff[, "Name"])
nm = setdiff(nm, "")

mat = matrix(0, nrow=length(nm), ncol=length(f))
rownames(mat)=nm
colnames(mat)=f

for(j in f)
  mat[ gff[ gff[,"feature"]==j, "Name"], j ] = 1

mat = mat[ rowSums(mat)>0, ]

z = apply(mat, 1, function(z) paste(f[as.logical(z)], collapse=":"))

sink("table-ncRNAs.txt")
cat(nrow(mat), "unique features ('uniqueness' defined via feature name).\n\n")
print(table(z))
sink()

pdf(file="table-ncRNAs.pdf", width=20, height=23)
heatmap(mat, col=c("white", "black"))
dev.off()

