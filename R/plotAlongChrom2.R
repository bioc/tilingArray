plotAlongChrom2 = function(chr, coord, highlight, segObj, y, ylim, probeAnno, isDirectHybe=FALSE,
  scoreShow="pt", nrBasesPerSeg, gff, haveNames=TRUE, haveLegend=TRUE, main="", colors,
  pointSize=unit(0.6, "mm")) {
  

  VP = c(title=0.2, expr1=5, z1=0.4, gff1=1, coord=1, gff2=1, z2=0.4, expr2=5, legend=0.4)

  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
    "cp" = "#101010", "highlight" = "red", "threshold" = "grey")
  if(!missing(colors)) {
    mt = match(names(colors), names(defaultColors))
    if(any(is.na(mt)))
      stop(paste("Cannot use color specification for", names(colors)[is.na(mt)]))
    defaultColors[mt] = colors 
  }
  colors = defaultColors
    
  ##if(!haveLegend)
  ##  VP = VP[-which(names(VP)=="legend")]

  ## do not draw p-value bars
  VP = VP[-which(names(VP)%in%c("z1", "z2"))]
  
  if(!missing(y)) {
    if(missing(probeAnno))
      stop("If 'y' is specified, 'probeAnno' must also be specified.")
    if(!missing(segObj))
      stop("If 'y' is specified, 'segObj' must not be specified.")
  } else {
    if(missing(segObj))
      stop("Please specify either 'y' or 'segObj'")
  }

  pushViewport(viewport(width=0.9, height=0.95)) ## plot margin
  pushViewport(viewport(layout=grid.layout(length(VP), 1, height=VP)))
  for(i in 1:2) {
    strand = c("+", "-")[i]
    threshold = as.numeric(NA)

    if(!missing(y)) {
      index = get(paste(chr, strand, "index", sep="."), envir=probeAnno)
      sta   = get(paste(chr, strand, "start", sep="."), envir=probeAnno)
      end   = get(paste(chr, strand, "end",   sep="."), envir=probeAnno)
      dat = list(mid = (sta+end)/2,
                 y   = y[index],
                 unique = get(paste(chr, strand, "unique", sep="."), envir=probeAnno))
      lengthChr = end[length(end)]
      sgs = NULL     
    } else {
      dat = get(paste(chr, strand, "dat", sep="."), segObj)
      dat$mid   = (dat[["start"]] + dat[["end"]])/2
      lengthChr = dat[["end"]][length(dat[["end"]])]
      
      if("threshold" %in% ls(segObj))
        threshold = get("threshold", segObj)
      
      if("segScore" %in% ls(segObj)) {
        sgs = get("segScore", segObj)
        if(!missing(nrBasesPerSeg))
          stop("Please do not specify 'nrBasesPerSeg' when 'segObj' contains 'segScore'")
      } else {
        if(missing(nrBasesPerSeg))
          stop("Please specify 'nrBasesPerSeg' ('segScore' was not found in 'segObj')")
        seg = get(paste(chr, strand, "seg", sep="."), segObj)
        cp  = round( lengthChr / nrBasesPerSeg)
        th  = c(1, seg$th[cp, 1:cp])
        sgs = data.frame(
          chr    = I(rep(chr, cp)),
          strand = I(rep(strand, cp)),
          start  = dat[["start"]][dat[["ss"]]][th[-length(th)]],
          end    = dat[["end"]][dat[["ss"]]][th[-1]]-1)
      }
    } 
    
    if(missing(coord))
      coord = c(1, lengthChr)
    if(missing(ylim))
      ylim = range(dat[["y"]])

   
    plotSegmentation(x= dat[["mid"]],
                     y= dat[["y"]], xlim=coord, ylim=ylim, uniq=dat$unique,
                     segScore=sgs, threshold=threshold, scoreShow=scoreShow,
                     gff=gff, chr=chr,
                     strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                     VP=VP, colors=colors, pointSize=pointSize, haveNames=haveNames)
    
  }

  ## chromosomal coordinates
  pushViewport(dataViewport(xData=coord, yscale=c(-0.4,0.8), extension=0, 
                            layout.pos.col=1, layout.pos.row=which(names(VP)=="coord")))
  grid.lines(coord, c(0,0), default.units = "native")
  tck = alongChromTicks(coord)
  grid.text(label=formatC(tck, format="d"), x = tck, y = 0.2, 
            just = c("centre", "bottom"), gp = gpar(cex=.6), default.units = "native")
  grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17,  default.units = "native")

  if(!missing(highlight)){
    mt = (match(highlight$strand, c("-", "+"))-1.5)*2
    co = highlight$coord
    if(is.na(mt) || !is.numeric(co))
      stop("Invalid parameter 'highlight'.")
    grid.segments(x0=co, x1=co, y0=c(0,0), y1=c(0.4,0.4)*mt, default.units = "native",
                  gp=gpar(col=colors["highlight"], lwd=2))
  }
  popViewport()

  ## title
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(VP)=="title")))
  grid.text(label=paste(main, "Chr", chr), x = 0.5, y = 0.5, 
            just = "centre", gp = gpar(cex=1))
  popViewport()

  ## legend
  if(haveLegend)
    plotAlongChromLegend(which(names(VP)=="legend"))
  
  popViewport(2)
}

## gff and chrSeqname into an environment or object?

plotSegmentation = function(x, y, xlim, ylim, uniq, segScore, threshold, scoreShow,
  gff, chr, strand, VP, colors, pointSize, haveNames, probeLength=25) {

  ## could this be done better?
  if(is.matrix(y))
    y = rowMeans(y)
    
  stopifnot(length(x)==length(y), length(x)==length(uniq))

  if(missing(xlim)) {
    xlim=range(x)
  } else {
    sel = (x>=xlim[1])&(x<=xlim[2])
    x = x[sel]
    y = y[sel]
    uniq = uniq[sel]
  }
  
  istrand = match(strand, c("+", "-"))
  stopifnot(length(strand)==1, !is.na(istrand))

  if(!is.na(threshold)) {
    y = y-threshold
    ylim = ylim-threshold
  }
  
  ## the expression data. use two viewports for different clipping behavior
  vpr = which(names(VP)==sprintf("expr%d", istrand))
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
    layout.pos.col=1, layout.pos.row=vpr))
  grid.yaxis()
  
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  ord  = c(which(!uniq), which(uniq))
  colo = ifelse(uniq[ord], colors[strand], colors["duplicated"])

  if(!is.na(threshold))
    grid.lines(y=unit(0, "native"), gp=gpar(col=colors["threshold"]))
  
  if(!is.null(segScore)) {
    segSel   = which(segScore$chr==chr & segScore$strand==strand)
    segstart = segScore[segSel, "start"]
    segend   = segScore[segSel, "end"]
    diffss =  segstart[-1] - segend[-length(segend)]
    meanss = (segstart[-1] + segend[-length(segend)])/2
    stopifnot(all(diffss>=(-probeLength)))
    
    grid.segments(x0 = unit(meanss, "native"), x1 = unit(meanss, "native"),
                  y0 = unit(0.1, "npc"),       y1 = unit(0.9, "npc"),
                  gp = gpar(col=colors["cp"]))

    grid.points(x[ord], y[ord], pch=20, size=pointSize, gp=gpar(col=colo))

  }
  popViewport(2)

  if(FALSE) {
  ## if(!is.null(segScore)) {
    pushViewport(dataViewport(xData=xlim, yscale=c(0,2), extension=0, clip="on",
       layout.pos.col=1, layout.pos.row=which(names(VP)==sprintf("z%d", istrand))))

    deckel = function(p, pmax=10) {
      p = -log(p, 10)
      p[p>pmax] = pmax
      sqrt(p/pmax)
    }
    
    if(scoreShow %in% colnames(segScore)) {
      colo  = rep("white", length(segSel))
      val   = deckel(segScore[segSel, scoreShow])
      isa   = !is.na(val)
      tmp   = colorRamp(brewer.pal(9, "Blues")[-9])(val[isa]) / 256
      colo[!is.na(val)] = rgb(tmp[,1], tmp[,2], tmp[,3])
    } else {
      colo = "#f0f0f0"
    }
    grid.rect(x = unit(segstart, "native"), 
              y = unit(0.4, "npc"),
              width  = unit(segend-segstart+1, "native"),
              height = unit(0.8, "npc"),
              just   = c("left", "center"),	
              gp     = gpar(col="#a0a0a0", fill=colo))
    popViewport(1)
  }
  
  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=which(names(VP)==sprintf("gff%d", istrand))))

  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  sel = which(gff[, "chr"] == chr &
              gff[, "strand"]  == strand &
              gff[, "start"] <= xlim[2] &
              gff[, "end"]   >= xlim[1])

  nam1    = gff[sel, "gene"]
  featnam = gff[sel, "Name"]
  featnam[!is.na(nam1)] = nam1[!is.na(nam1)]
  featsp  = split(seq(along=sel), gff$feature[sel])

  ### gene ###
  whnames = integer(0)
  if("gene" %in% names(featsp)) {
    i = featsp[["gene"]]
    s = sel[i]
    grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, y1 = 0,
                  default.units = "native", gp = gpar(col="#a0a0a0"))
    whnames = i
    ## cat(paste(featnam[i], gff$start[s], gff$end[s], sep="\t", collapse="\n"), "\n\n")
  }
  
  featDraw = featureDrawing()
  sfeatsp  = featsp[rownames(featDraw)]
  ll       = listLen(sfeatsp)

  if(any(ll>0)) {
    i      = unlist(sfeatsp)
    ord    = order(gff$start[sel[i]])
    gp     = gpar(col = rep(featDraw$col,  ll)[ord],
                 fill = rep(featDraw$fill, ll)[ord])
    i      = i[ord]
    s      = sel[i]
    
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[ names(sfeatsp)!="CDS" ]))
  }

  if(haveNames && (length(whnames)>0)) {
    txtcex = 0.7
    txtdy  = 0.7
    whnames = whnames[!duplicated(featnam[whnames])]
    s      = sel[whnames]
    txtx   = (gff$start[s]+gff$end[s])/2
    txty   = numeric(length(s))
    ord    = order(txtx)
    whnames = whnames[ord]
    s      = s[ord]
    txtx   = txtx[ord]
    
    strw   = convertWidth(stringWidth(featnam[whnames]), "native", valueOnly=TRUE)*txtcex
    rightB = txtx[1] + 0.5*strw[1]
    doText = rep(TRUE, length(whnames))
    if(length(whnames)>1) {
      for(k in 2:length(whnames)) {
        leftB = txtx[k] - 0.5*strw[k]
        if(leftB > rightB) {
          rightB = txtx[k] + 0.5*strw[k]
        } else {
          if(!any(txty[k-(1:2)]==txtdy)) {
            txty[k]= txtdy
          } else {
            if(!any(txty[k-(1:2)]== -txtdy)) {
              txty[k]= -txtdy
            } else {
              doText[k] = FALSE
            }
          }
        } ##  else
      } ## for
    }
    grid.text(label = featnam[whnames][doText],
              x = txtx[doText], y = txty, gp=gpar(cex=txtcex), 
              default.units = "native")
  } ## if
  
  ### intron ###
  if("intron" %in% names(featsp)) {
    i = featsp[["intron"]]
    s = sel[i]
    mid = (gff$start[s]+gff$end[s])/2
    wid = (gff$end[s]-gff$start[s])/2 
    for(z in c(-1,1))
      grid.segments(x0 = mid,
                    x1 = mid+z*wid,
                    y0 = c(1, -1)[istrand]*1.2,  ## istrand is 1 or 2
                    y1 = c(1, -1)[istrand]*0.95,
                    default.units = "native",
                    gp = gpar(col="black"))
    ## cat(paste(featnam[i], gff$start[s], gff$end[s], sep="\t", collapse="\n"), "\n\n")
  } ## if
  

  ##notDealtWith = setdiff(names(featsp), c(rownames(featDraw), "gene", "intron"))
  ##if(length(notDealtWith)>0)
  ##  cat("Not displayed:", paste(notDealtWith, collapse=", "), "\n")
   
  popViewport()
}

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks = function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/5, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}

##------------------------------------------------------------
## legend
##------------------------------------------------------------
plotAlongChromLegend = function(vpr) {
  featDraw = featureDrawing()

  featDraw = featDraw[!(rownames(featDraw) %in% c("transposable_element_gene", "transposable_element")), ]
  dx       = 1/(1+nrow(featDraw))
  i        = 1:nrow(featDraw)

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr))

  grid.lines(c(0,1), c(1,1), default.units = "npc")

  
  grid.rect(x     = (i-1)*dx,
            y     = 0.4,
            width = dx*0.3,
            height= 0.8, 
            default.units = "npc", just  = c("left", "center"),
            gp    = do.call("gpar", featDraw))

  grid.text(label = rownames(featDraw),
            x     = (i-0.65)*dx,
            y     = 0.4,
            default.units = "npc", just  = c("left", "center"),
            gp    = gpar(cex=.7))

  popViewport()
}


##------------------------------------------------------------
## plotDuplication
##------------------------------------------------------------
plotDuplication = function(xlim, chr, strand, probeAnno, VP) { 

  istrand = match(strand, c("+", "-"))
  stopifnot(length(strand)==1, !is.na(istrand))

  sta = get(paste(chr, strand, "start", sep="."),  probeAnno)
  uni = get(paste(chr, strand, "unique", sep="."), probeAnno)
  ord = order(sta)
  uni = uni[ord]
  sta = sta[ord]

  pp  = diff(uni)
  x0 = which(pp==-1)  ## transition from TRUE to FALSE
  x1 = which(pp==+1)  ## transition from FALSE to TRUE

  if(!uni[1])
    x0 = c(1, x0)
  if(!uni[length(uni)])
    x1 = c(x1, length(uni))

  stopifnot(length(x0)==length(x1))
  stopifnot(all(x1>x0))
  
  pushViewport(dataViewport(xData=xlim, yscale=c(-1,1),  extension=0,  clip="on",
    layout.pos.col=1, layout.pos.row=which(names(VP)==sprintf("dup%d", istrand), "id")))

  grid.segments(x0=sta[x0], x1=sta[x1], y0=0, y1=0, default.units = "native",
                gp=gpar(lwd=3, col="#606060")) 
  popViewport()
}


##------------------------------------------------------------
## featureDrawing
##------------------------------------------------------------
featureDrawing = function() {
  res = data.frame(
    col      = I(c("#d0e0f0", "#d94801", "#005a32", "#707070", "#fc4e2a", 
                   "#A65628", "#A65628", "#FED976", "#FED976")),
    fill     = I(c("#7AADD1", "#fd8d3c", "#41ab5d", "#e0e0e0", "#feb24c", 
                   "#BF5B17", "#BF5B17", "#FFEDA0", "#FFEDA0"))) 
  rownames(res) =   c("CDS",     "tRNA",    "snoRNA",  "pseudogene", "ncRNA",   
            "transposable_element", "transposable_element_gene", "centromere", "telomere")
  return(res)
}
