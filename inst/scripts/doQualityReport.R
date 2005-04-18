library("tilingArray")
library("affy")

options(error=recover)

if(!exists("probeAnno"))
  load("probeAnno.rda")

normRefFiles = file.path("Celfiles",
  c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz",
    "041119_S96genDNA_re-hybe.cel.gz",
    "041120_S96genDNA_re-hybe.cel.gz"))

files = file.path("Celfiles",
  c("050415_totcDNA_20ug_Affy11.cel.gz",
    "050415_totcDNA_20ug_Prok16.cel.gz"))

hybeType=c(rep("Reverse", 2))  ## "Direct"

if(!exists("x"))
  x = read.affybatch(filenames=files, compress=TRUE, verbose=TRUE)

if(!exists("normRef"))
  normRef = read.affybatch(filenames=normRefFiles, compress=TRUE, verbose=TRUE)

for(i in seq(along=files))
  qualityReport(x = x[,i], hybeType=hybeType[i], normRef = normRef,
              selectGenes = c("YLR342W", "YOR153W", "YDL185W", "YKL182W",
                "YGR027C", "YDR384C"),
              gff = gff, probeAnno = probeAnno,
              outputDir = "qualityReports")

