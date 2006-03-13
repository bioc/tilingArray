plotAlongChrom = function(segObj, y, probeAnno, gff,
                          nrBasesPerSeg,
                          isDirectHybe=FALSE, 
                          what = c("dots"), ## "heatmap"
                          chr, coord, highlight, ylim, 
                          colors, 
                          doLegend=TRUE,
                          featureColorScheme=1,
                          featureExclude=c("chromosome","gene","nucleotide_match", "insertion", "intron"),
                          featureNoLabel=c("uORF"),
                          pointSize=unit(0.6, "mm"),
                          main, ...) {

  ## set up the viewports of the plot layout.
  VP = c("title"=0.8, "expr+"=5, "gff+"=1, "coord"=1, "gff-"=1, "expr-"=5, "legend"=0.4)
  if(!doLegend)
     VP = VP[-which(names(VP)=="legend")]

  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
    "cp" = "#777777", "highlight" = "red", "threshold" = "grey")
  if(!missing(colors)) {
    mt = match(names(colors), names(defaultColors))
    if(any(is.na(mt)))
      stop(paste("Cannot use color specification for", names(colors)[is.na(mt)]))
    defaultColors[mt] = colors 
  }
  colors = defaultColors

  ## check that either y or segObj is present
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

    ## extract and treat  the data
    threshold = as.numeric(NA)

    ## three mutually exclusive cases:
    ## 1.) y and probeAnno
    ## 2.) segObj is an environment
    ## 3.) segObj is object of S4 class "segmentation"
    if(!missing(y)) {
      index = get(paste(chr, strand, "index", sep="."), envir=probeAnno)
      sta   = get(paste(chr, strand, "start", sep="."), envir=probeAnno)
      end   = get(paste(chr, strand, "end",   sep="."), envir=probeAnno)
      dat = list(x = (sta+end)/2,
                 y = y[index,, drop=FALSE],
                 flag = get(paste(chr, strand, "unique", sep="."), envir=probeAnno))
      stopifnot(is.numeric(dat$flag))
      lengthChr = end[length(end)]
      
    } else {
      if(inherits(segObj, "segmentation")){
        ## new: S4 class
        dat$x = segObj@x
        dat$y = segObj@y
        dat$flag = 123456
        browser()
      } else {
        ## old: environment
        dat = get(paste(chr, strand, "dat", sep="."), segObj)
        stopifnot(all(c("start", "end", "unique", "ss") %in% names(dat)))
        dat$x = (dat$start + dat$end)/2
        dat$flag = dat$unique
        lengthChr = dat$end[length(dat$end)]
        
        if("theThreshold" %in% ls(segObj))
          threshold = get("theThreshold", segObj)
        
        if("segScore" %in% ls(segObj)) {
          if(!missing(nrBasesPerSeg))
            stop("Please do not specify 'nrBasesPerSeg' when 'segObj' contains 'segScore'")
          sgs = get("segScore", segObj)
          
        } else {
          if(missing(nrBasesPerSeg))
            stop("Please specify 'nrBasesPerSeg' ('segScore' was not found in 'segObj')")
          seg = get(paste(chr, strand, "seg", sep="."), segObj)
          cp  = round( lengthChr / nrBasesPerSeg)
          th  = c(1, seg$th[cp, 1:cp])
          sgs = list(start  = dat$start[dat$ss][th[-length(th)]],
                     end    = dat$end[dat$ss][th[-1]]-1)
        }
        dat$estimate = (sgs$start[-1] + sgs$end[-length(sgs$end)]) / 2
      } 
    }
    ## At this point, no matter what the input to the function was,
    ##    we have the list 'dat' with elements
    ## x: x coordinate
    ## y: y coordinate
    ## flag
    ## estimate (optional)
    ## lower (optional)
    ## upper (optional)
    
    if(missing(coord))
      coord = c(1, lengthChr)
    
    ## plot the data
    vpr=which(names(VP)==sprintf("expr%s", strand))
    switch(what,
      "dots" = {
        if(missing(ylim))
          ylim = quantile(dat$y[dat$x>=coord[1] & dat$x<=coord[2]], c(0.01, 0.99), na.rm=TRUE)
   
        plotSegmentationDots(dat, xlim=coord, ylim=ylim, 
                             threshold=threshold, 
                             chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                             vpr=vpr, colors=colors, pointSize=pointSize, ...)
      },
      "heatmap" = {
        plotSegmentationHeatmap(dat, xlim=coord, 
                                threshold=threshold, 
                                chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                                vpr=vpr, colors=colors, ...)
      },
           stop(sprintf("Invalid value '%s' for argument 'what'", what))
    ) ## switch
  
    ## plot the features
    plotFeatures(gff=gff, chr=chr, xlim=coord, strand=strand,
                 vpr=which(names(VP)==sprintf("gff%s", strand)), 
                 featureColorScheme=featureColorScheme,
                 featureExclude=featureExclude, featureNoLabel=featureNoLabel)
  
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
  grid.text(label=paste("Chr ", chr, sep=""), x=0.5, y=1, just="centre", gp=gpar(cex=1))
  if(!missing(main))
    grid.text(label=main, x=-0.1, y=1, just="left", gp=gpar(cex=1))
  popViewport()

  ## legend
  if(doLegend)
    plotAlongChromLegend(which(names(VP)=="legend"),
         featureColorScheme=featureColorScheme, featureExclude=featureExclude)
  
  popViewport(2)
}

## ------------------------------------------------------------
## plot Segmentation with Dots
## ------------------------------------------------------------
plotSegmentationDots = function(dat, xlim, ylim, threshold, 
  chr, strand, vpr, colors, pointSize) {

  if(is.matrix(dat$y))
    dat$y = rowMeans(dat$y) ##  if >1 samples, take mean over samples
  stopifnot(length(dat$y)==length(dat$x), length(dat$flag)==length(dat$x))
  
  if(missing(xlim)) {
    xlim=range(dat$x, na.rm=TRUE)
  } else {
    sel = (dat$x>=xlim[1])&(dat$x<=xlim[2])
    dat$x = dat$x[sel]
    dat$y = dat$y[sel]
    dat$flag = dat$flag[sel]
  }
  
  if(!is.na(threshold)) {
    dat$y = dat$y-threshold
    ylim = ylim-threshold
  }
  
  ## the expression data. use two viewports for different clipping behavior
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
    layout.pos.col=1, layout.pos.row=vpr))
  grid.yaxis()
  
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  ord  = c(which(dat$flag!=0), which(dat$flag==0))
  colo = ifelse(dat$flag[ord]==0, colors[strand], colors["duplicated"])

  if(!is.na(threshold))
    grid.lines(y=unit(0, "native"), gp=gpar(col=colors["threshold"]))

  ## segment boundaries
  if(!is.null(dat$estimate)) {
    grid.segments(x0 = unit(dat$estimate, "native"),
                  x1 = unit(dat$estimate, "native"),
                  y0 = unit(0.1, "npc"),
                  y1 = unit(0.9, "npc"),
                  gp = gpar(col=colors["cp"]))
  }
  
  grid.points(dat$x[ord], dat$y[ord], pch=20, size=pointSize, gp=gpar(col=colo))
  popViewport(2)

} ## plotSegmentationDots

##------------------------------------------------------------- 
##  plotSegmentationHeatmap
##------------------------------------------------------------- 
plotSegmentationHeatmap = function(dat, xlim, threshold,
  chr, strand, vpr, colors, transformation=function(z) z) {

  if(missing(xlim)) {
    xlim=range(dat$x, na.rm=TRUE)
  } else {
    sel = (dat$x>=xlim[1])&(dat$x<=xlim[2])
    dat$x = dat$x[sel]
    dat$y = dat$y[sel,, drop=FALSE ]
    dat$flag = dat$flag[sel]
  }

  ord = order(dat$x)
  dat$x = dat$x[ord]   ## sort by x-coordinates to simplify smoothing
  dat$y = dat$y[ord,, drop=FALSE]
  dat$flag = dat$flag[ord]
  
  ## Use two viewports for different clipping behavior
  ylim = c(-1, 2+ncol(dat$y))
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
    layout.pos.col=1, layout.pos.row=vpr))

  ylab = colnames(dat$y)
  grid.yaxis( (1:ncol(dat$y)), ylab, gp=gpar(cex=0.5))

  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  ord  = c(which(dat$flag!=0), which(dat$flag==0))
  colo = ifelse(dat$flag[ord]==0, colors[strand], colors["duplicated"])

  grid.image(dat$x, 1:ncol(dat$y), z=matrix(transformation(dat$y), ncol=ncol(dat$y), nrow=nrow(dat$y)),
             xlim=xlim, uniq=dat$flag)

  ## segment boundaries
  if(!is.null(dat$estimate)) {
    grid.segments(x0 = unit(dat$estimate, "native"),
                  x1 = unit(dat$estimate, "native"),
                  y0 = unit(0.1, "npc"),
                  y1 = unit(0.9, "npc"),
                  gp = gpar(col=colors["cp"]))
  }

  popViewport(2)
} ## end of plotSegmentationHeatmap


## ------------------------------------------------------------
## plot Features
## ------------------------------------------------------------
plotFeatures = function(gff, chr, xlim, strand, vpr, featureColorScheme, featureExclude, featureNoLabel) {

  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  sel = which(gff[, "chr"] == chr &
              gff[, "strand"]  == strand &
              gff[, "start"] <= xlim[2] &
              gff[, "end"]   >= xlim[1])

  geneName = gff[sel, "gene"]
  featName = gff[sel, "Name"]
  featName[!is.na(geneName)] = geneName[!is.na(geneName)]

  feature  = as.character(gff[sel, "feature"])

  ## split by feature type (e.g. CDS, ncRNA)
  featsp = split(seq(along=sel), feature)

  ## drop the ignorable ones
  featsp = featsp[ ! (names(featsp) %in% featureExclude) ]
  
  ## gene: just a horizontal line
  whnames = integer(0)
  if("gene" %in% names(featsp)) {
    i = featsp[["gene"]]
    s = sel[i]
    grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, y1 = 0,
                  default.units = "native", gp = gpar(col="#a0a0a0"))
    whnames = i
  }

  ## colors for boxes
  featCols = featureColors(featureColorScheme)

  ## check that we know how deal with all features
  whm = names(featsp) %in% rownames(featCols)
  if(!all(whm))
    stop("Don't know how to handle feature(s) '", paste(names(featsp)[!whm], collapse=", "), "'.", sep="")

  sfeatsp  = featsp[rownames(featCols)]
  ll       = listLen(sfeatsp)

  if(any(ll>0)) {
    i  = unlist(sfeatsp)
    gp = gpar(col = rep(featCols$col,  ll),
                 fill = rep(featCols$fill, ll))
    s  = sel[i]
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% featureNoLabel)]))
    ## additional potentially useful values for featureNoLabel: "binding_site", "TF_binding_site"
  }

  ## labels
  if( !all(tolower(featureNoLabel)=="all") && (length(whnames)>0)) {

    ## this is a bit of a hack to abbreviate the labels of "binding site" features:
    bindingRegexpr = "binding.?site.*$"
    isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 0)
    if(any(isBindingSite)) {
      ## replace long labels
      featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
    }

    ## remove duplicated names that are not binding sites
    whnames = whnames[isBindingSite | !duplicated(featName[whnames])]

    txtcex = 0.6
    txtdy  = 0.7
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
    if(length(whnames) >1) {
      for(k in 2:length(whnames)) {
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
  } ## if
  popViewport()

} ## plotFeatures

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks = function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
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
plotAlongChromLegend = function(vpr, nr=2, 
    featureColorScheme, featureExclude){
  
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
  

  featCols = featureColors(featureColorScheme)
  featCols = featCols[ !(rownames(featCols) %in% featureExclude), ]

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr, yscale=c(0.5, nr+0.5)))

  i = 1:nrow(featCols)
  for(r in 1:nr)
    formatRow(featCols[ceiling(i/nrow(featCols)*nr-1e-10)==r, ], row=nr-r+1)
  
  popViewport()
  
} ## plotAlongChromLegend

##------------------------------------------------------------
## featureColors
## note that features are drawn in the order in which they appear
## here, this can be used to let important features overdraw less
## important ones (e.g. tRNA is more specific than ncRNA)
## to test, say tilingArray:::plotAlongChromLegend()
##------------------------------------------------------------
featureColors = function(scheme=1){
  
  defaultColors = c(
    "chromosome"  = NA,
    "nucleotide_match" = "#e0e0e0",   ## light gray
    "pseudogene"  = "#e0e0e0",        ## light gray
    "uORF"        =   "#FED976" ,     ## orange
    "nc_primary_transcript" = "#a0a0a0",    ## grey
    "region" = "#cc66cc",           ## light red-violet	
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
    "ncRNA"       = "#a6d96a",    ## green 
    "tRNA"        = "#a6d96a",    ## green
    "snRNA"       = "#8C6BB1",    ## purple
    "rRNA"        = "#fdae61",    ## meat
    "snoRNA"      = "#7F5A58",    ## red brown
    "binding_site"    = "#C9C299", ## lemon chiffon
    "TF_binding_site" = "#C9C299" ## lemon chiffon
  )
  darkenborder = as.logical(c(rep(1,3),0,rep(1, 17),0,0))
  stopifnot(length(darkenborder)==length(defaultColors))
  
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
  border = ifelse(darkenborder, darken(fill), fill)
  
  res = data.frame(fill=I(fill),
    col =I(border))
  rownames(res)=names(defaultColors) 
  return(res)
}
