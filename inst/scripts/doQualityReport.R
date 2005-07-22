library("tilingArray")
library("affy")

options(error=recover)

if(!exists("probeAnno"))
  load("probeAnno.rda")
if(!exists("a"))
  load("a.rda")

hybeSets = list(
  "polyA2" = c("05_04_27_2xpolyA_NAP3.cel.gz",
    "05_04_26_2xpolyA_NAP2.cel.gz",
    "05_04_20_2xpolyA_NAP_2to1.cel.gz"),
  "tot" = c("050409_totcDNA_14ug_no52.cel.gz",
    "030505_totcDNA_15ug_affy.cel.gz"),
  "dir" = c("050621_dirPolyARNA_10ug_2-3.cel.gz",
    "050621_dirPolyARNA_10ug_2-3_4x.cel.gz"))

hsAnno = data.frame(
  dir = I(c("seg-polyA-050525", "seg-tot-050525", "seg-dir-050721")),
  isDirect = c(FALSE, FALSE, TRUE))
rownames(hsAnno) = c("polyA2", "tot", "dir")

outdir = "qualityReports"

##--------------------------------------------------
makeFig1 = function(y, main="") {

  ## push 1: first level of layout
  pushViewport(viewport(layout=grid.layout(2, 1, height=c(1.5, 2), width=1)))

  ##
  ## A) 14:380-480: overview
  ##
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,
                        layout=grid.layout(1, 2, height=1, width=c(0.01, 1))))
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
  plotAlongChrom(chr=14, coord= c(410000,488150), y=y,
                  colors=c(cp="#d0d0d0"), main=main,
                  probeAnno = probeAnno, gff=gff,
                  haveNames=FALSE, haveLegend=FALSE, pointSize=unit(0.1, "mm"),
                  featColScheme=2)
  popViewport(2)

  ## second level of layout
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=2,
                        layout=grid.layout(2, 6, height=c(1, 1), width=c(0.1, 1, 0.07, 1, 0.07, 1))))
  
  ##
  ## B) 13:550k splicing RPS16A, RPL13B
  ## middle row, left
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
  plotAlongChrom(chr=13, coord = c(550044, 553360), y=y, 
                  probeAnno = probeAnno, gff=gff, haveLegend=FALSE, colors=c(cp="#d0d0d0"))
  popViewport()

  ## C) 11:65k: MNN4, novel architecture
  ## middle row, middle
  pushViewport(viewport(layout.pos.col=4, layout.pos.row=1))
  plotAlongChrom(chr=11, coord = c(63370, 68270), y=y, 
                  probeAnno = probeAnno, gff=gff, haveLegend=FALSE, colors=c(cp="#d0d0d0"))
  popViewport()
  
  ## D) overlapping transcripts
  ## middle row, right
  pushViewport(viewport(layout.pos.col=6, layout.pos.row=1))
  plotAlongChrom(chr=14, coord = c(342500, 347545), y=y, 
                  probeAnno = probeAnno, gff=gff, haveLegend=FALSE, colors=c(cp="#d0d0d0"))
  popViewport()
  
  ## E) 2:360.5-366.5: novel isolated
  ## bottom row, left
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
  plotAlongChrom(chr=2, coord = c(360500, 365970), y=y, 
                  probeAnno = probeAnno, gff=gff, haveLegend=FALSE, colors=c(cp="#d0d0d0"))
  popViewport()

  ## F) 9:221-227: novel antisense SPO22
  ## bottom row, middle
  pushViewport(viewport(layout.pos.col=4, layout.pos.row=2))
  plotAlongChrom(chr=9, coord = c(221000, 226500), y=y, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE, colors=c(cp="#d0d0d0"))
  popViewport()
  
  ## pop 2,1
  popViewport(2)
  
}
## --------------------------------------------------

for(hs in seq(along=hybeSets)) {
  files = hybeSets[[hs]]
  load(file.path(hsAnno$dir[hs], "xn.rda"))
  for(fn in files) {
    f1 = file.path(outdir, sub(".cel.gz", ".pdf", fn))
    f2 = file.path(outdir, sub(".cel.gz", ".tiff", fn))
    pdf(file=f1, width=8, height=12)
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2, 1, height=c(1,1), width=1)))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    makeFig1(log(exprs(a)[,fn],2), main=paste("raw (log2 scale)", sub(".cel.gz", "", fn), "- "))
    popViewport(1)
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
    makeFig1(exprs(xn)[,fn], main=paste("normalized", sub(".cel.gz", "", fn), "- "))
    popViewport(2)
    dev.off()

    cmd = paste("convert -density 180", f1, "-compress RLE", f2, "&")
    cat(cmd, "\n")
    system(cmd)
  }
}
