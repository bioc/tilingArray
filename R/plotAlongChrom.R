plotAlongChrom = function(chr, coord, highlight, segObj, y, ylim, nrBasesPerSeg, 
                      probeAnno, gff,
                      colors, featColScheme=1,
                      isDirectHybe=FALSE, scoreShow = "pt", 
                      haveNames=TRUE, haveLegend=TRUE, main="", 
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
      stopifnot(is.numeric(dat$unique))
      lengthChr = end[length(end)]
      sgs = NULL     
    } else {
      dat = get(paste(chr, strand, "dat", sep="."), segObj)
      dat$mid   = (dat[["start"]] + dat[["end"]])/2
      lengthChr = dat[["end"]][length(dat[["end"]])]
      
      if("theThreshold" %in% ls(segObj))
        threshold = get("theThreshold", segObj)

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

    px = dat[["mid"]]
    py = dat[["y"]]
    
    if(missing(ylim))
      ylim = quantile(py[px>=coord[1] & px<=coord[2]], c(0.01, 0.99), na.rm=TRUE)
   
    plotSegmentation(x=px, y=py, xlim=coord, ylim=ylim, uniq=dat$unique,
                     segScore=sgs, threshold=threshold, scoreShow=scoreShow,
                     gff=gff, chr=chr,
                     strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                     VP=VP, colors=colors, pointSize=pointSize, haveNames=haveNames,
                     featColScheme=featColScheme)
    
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
  grid.text(label=paste(main, "Chr ", chr, sep=""), x = 0.5, y = 0.5, 
            just = "centre", gp = gpar(cex=1))
  popViewport()

  ## legend
  if(haveLegend)
    plotAlongChromLegend(which(names(VP)=="legend"))
  
  popViewport(2)
}

## gff and chrSeqname into an environment or object?

plotSegmentation = function(x, y, xlim, ylim, uniq, segScore, threshold, scoreShow,
  gff, chr, strand, VP, colors, pointSize, haveNames, probeLength=25, featColScheme,
  noTypeLabel = c("CDS"), #,"binding_site", "TF_binding_site"),
  exclude=c("chromosome","gene","nucleotide_match","insertion","intron")
    #   ,"ARS","repeat_region","repeat_family",  # new version: 2005-08-30 J
) {
  ## could this be done better?
  if(is.matrix(y))
    y = rowMeans(y) #  if >1 samples? -> take mean over samples
    
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

  ord  = c(which(uniq!=0), which(uniq==0))
  colo = ifelse(uniq[ord]==0, colors[strand], colors["duplicated"])

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
  }
  grid.points(x[ord], y[ord], pch=20, size=pointSize, gp=gpar(col=colo))
  popViewport(2)


  ## This code would plot a false color representation of some segment score
  ## (e.g. "p-value")
  if(FALSE) {
    if(!is.null(segScore)) {
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
  }

  ## features
  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=which(names(VP)==sprintf("gff%d", istrand))))

  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  sel = which(gff[, "chr"] == chr &
              gff[, "strand"]  == strand &
              gff[, "start"] <= xlim[2] &
              gff[, "end"]   >= xlim[1])

  geneName = gff[sel, "gene"]
  featName = gff[sel, "Name"]
  featName[!is.na(geneName)] = geneName[!is.na(geneName)]

  feature  = as.character(gff[sel, "feature"])

  ## split!
  featsp  = split(seq(along=sel), feature)

  ### gene ###
  whnames = integer(0)
  if("gene" %in% names(featsp)) {
    i = featsp[["gene"]]
    s = sel[i]
    grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, y1 = 0,
                  default.units = "native", gp = gpar(col="#a0a0a0"))
    whnames = i
  }

  featCols = featureColors(featColScheme, exclude)

  ## check that we know how deal with all features
  whm = names(featsp) %in% c(rownames(featCols), exclude) # change: 2005-08-25 J
  if(!all(whm))
    stop("Don't know how to handle feature(s) '", paste(names(featsp)[!whm], collapse=", "), "'.", sep="")

  sfeatsp  = featsp[rownames(featCols)]
  ll       = listLen(sfeatsp)

  if(any(ll>0)) {
    i      = unlist(sfeatsp)
    # change, 05/09/10 J, do ordering  by specificity not by start base
    # ord    = order(gff$start[sel[i]])
    ord = 1:length(i)

    gp     = gpar(col = rep(featCols$col,  ll)[ord],
                 fill = rep(featCols$fill, ll)[ord])
    i      = i[ord]
    s      = sel[i]
    
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% noTypeLabel)]))
    #if (length(grep("binding", names(sfeatsp)))>0) browser()
  }

  if(haveNames && (length(whnames)>0)) {
    bindingRegexpr = "binding.?site.*$"
    isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 0)
    if(any(isBindingSite)) {
      ## replace long labels
      featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
    }
    ## remove duplicated names that are not binding sites
    whnames = whnames[isBindingSite | !duplicated(featName[whnames])]

    txtcex = 0.7
    txtdy  = 0.7
    nNames  <- length(whnames) 
    s      = sel[whnames]
    txtx   = (gff$start[s]+gff$end[s])/2
    txty   = numeric(length(s))
    ord    = order(txtx)
    whnames = whnames[ord]
    s      = s[ord]
    txtx   = txtx[ord]

    
    strw   = convertWidth(stringWidth(featName[whnames]), "native", valueOnly=TRUE)*txtcex
    rightB = txtx[1] + 0.5*strw[1]
    doText = rep(TRUE, length(whnames))

    # not used so far:
    # textstarts <-  txtx - 0.5*strw
    # textends   <-  txtx + 0.5*strw

    # adjust text labels to be still readable in feature-dense areas:
    if(nNames>1) {
      for(k in 2:nNames) {
        leftB = txtx[k] - 0.5*strw[k]
        if(leftB > rightB) { # all texts not overlapping next to each other?
          rightB = txtx[k] + 0.5*strw[k]
        } else { # any overlaps?
          if(!any(txty[k-(1:2)]==txtdy)) {#  2 previous labels not moved up?
            txty[k]= txtdy                #   then this one 
          } else {                        #  else try move down:
            if(!any(txty[k-(1:2)]== -txtdy)) { 
              txty[k]= -txtdy             #  if 2 previous ones weren't
            } else {
              doText[k] = FALSE           #  otherwise don't put the label
            }
          }
        } ##  else
      } ## for
    }
    
    grid.text(label = featName[whnames][doText],
              x = txtx[doText], y = txty[doText], gp=gpar(cex=txtcex), 
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
    ## cat(paste(featName[i], gff$start[s], gff$end[s], sep="\t", collapse="\n"), "\n\n")
  } ## if
  

  ##notDealtWith = setdiff(names(featsp), c(rownames(featCols), "gene", "intron"))
  ##if(length(notDealtWith)>0)
  ##  cat("Not displayed:", paste(notDealtWith, collapse=", "), "\n")
   
  popViewport()
}

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks = function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/4, 10)
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
plotAlongChromLegend = function(vpr=1, nr=2, # was nr=3
    exclude=c("chromosome","gene","nucleotide_match","insertion","intron")
    #"repeat_region","repeat_family","ARS","intron","nc_primary_transcript") # change Aug 30,2005 J
  ) {
  formatRow = function(featColsOneRow, row) {
    ## print(featColsOneRow)
    strWid   = convertWidth(stringWidth(rownames(featColsOneRow)), "npc", valueOnly=TRUE)
    n        = length(strWid)
    inbetWid = 0.2*min(strWid)
    totWid   = sum(strWid)+(n-1)*inbetWid
    x        = c(0, cumsum(strWid[-n])) + (0:(n-1))*inbetWid 
    y        = numeric(length(x))

    x      = x/totWid
    strWid = strWid/totWid
    grid.rect(x = x, width = strWid, 
              y = unit(row, "native"), height = unit(1, "native")- unit(1, "mm"), 
              just  = c("left", "center"), default.units="npc",
              gp    = do.call("gpar", featColsOneRow))
    
    grid.text(label = rownames(featColsOneRow),
              x = unit(x + strWid/2, "native"), y = unit(row, "native"),
              just  = c("center", "center"), gp=gpar(cex=0.66))
  }
  

  featCols = featureColors(1,exclude)

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr, yscale=c(0.5, nr+0.5)))
  ## grid.lines(c(0,1), c(1,1), default.units = "npc", gp=gpar(col="black",lty=2))

  i = 1:nrow(featCols)
  for(r in 1:nr)
    formatRow(featCols[ceiling(i/nrow(featCols)*nr-1e-10)==r, ], row=nr-r+1)
  
  popViewport()
  
}#plotAlongChromLegend

##------------------------------------------------------------
## featureColors
## note that features are drawn in the order in which they appear
## here, this can be used to let important features overdraw less
## important ones (e.g. tRNA is more specific than ncRNA)
## to test, say tilingArray:::plotAlongChromLegend()
##------------------------------------------------------------
featureColors = function(scheme=1, exclude=c()) {

  defaultColors = c("chromosome"  = NA,
                    "nucleotide_match" = "#e0e0e0",    ## light gray
                    "pseudogene"  = "#e0e0e0",    ## light gray
                    "uORF"        = "#e0e0e0",    ## light gray
                    "nc_primary_transcript" = "#a0a0a0",    ## grey
			  "region" = "#cc66cc",    ## light red-violet	
                    "repeat_family" = "#CC6666",    ## light red
                    "repeat_region" = "#e31a1c",    ## bright red                    
                    "transposable_element"  = "#f1b6da",    ## pink
                    "transposable_element_gene"= "#f1b6da",   				      
                    "ARS"         = "#CC9966",    ## light brown
                    "centromere"  = "#FFEDA0",    ## orange
                    "telomere"    = "#FFEDA0",    ## orange
                    "insertion"   = "#FFEDA0",    ## orange
                    "CDS"         = "#addfff",    ## light blue
                    "CDS_dubious" = "#e0f1f2",    ## lighter blue
                    "ncRNA"       = "#3b9c9c",    ## cyan
                    "tRNA"        = "#a6d96a",    ## green
                    "snRNA"       = "#8C6BB1",    ## purple
                    "rRNA"        = "#fdae61",    ## meat
                    "snoRNA"      = "#7F5A58",    ## red brown
                    "binding_site"    = "#C9C299", ## lemon chiffon
                    "TF_binding_site" = "#C9C299") ## lemon chiffon

  # kick out unwanted Features, no need to return a color defintion for them
  defaultFeatures  <- names(defaultColors)
  selectedFeatures <- setdiff(defaultFeatures, exclude)
  defaultColors    <- defaultColors[selectedFeatures]

  fill = switch(scheme,
    default  = defaultColors,
    unicolor = ifelse(is.na(defaultColors), NA,  "#addfff"),  ## light blue
    stop("Sapperlot"))

  ## calculate hex string for a color that is a little bit darker than the
  ## hex string in the argument
  darken = function(x) {
    sel = !is.na(x)
    xRGB = hex2RGB(x[sel])
    xRGB@coords = 0.5 * coords(xRGB) ## unfortunately there is no coords<- method
    res = rep(as.character(NA), length(x))
    res[sel] = hex(xRGB)
    return(res)
  }
    
  res = data.frame(fill=I(fill),
                   col =I(darken(fill)))
  rownames(res)=names(defaultColors) 
  return(res)
}

