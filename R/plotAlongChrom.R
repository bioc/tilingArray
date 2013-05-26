# plotting along chromosome (extra bar/s, gff optional, legend at bottom optional)
# thresh argument missing in plotAlongChrom (used for setting global ylim for heatmaps)
# also need ability to change the color scheme for panel and main heatmap
plotAlongChrom = function(segObj, y, probeAnno, gff,
                          isDirectHybe=FALSE, 
                          what = c("dots"), ## "heatmap"
                          chr, coord, highlight,  
                          colors, 
                          doLegend=FALSE,
                          featureExclude=c("chromosome", "nucleotide_match", "insertion"),
                          featureColorScheme=1, extras, 
                          rowNamesHeatmap, rowNamesExtras, ylab, ylabExtras, main,
                          colHeatmap=colorRamp(brewer.pal(9, "YlGnBu")),
                          colExtras=colorRamp(brewer.pal(9, "Reds")), 
                          sepPlots=FALSE, reOrder=TRUE,...) {

  ## set up the viewports of the plot layout.
  VP = c("title"=0.4, "expr+"=5, "gff+"=1, "coord"=1, "gff-"=1, "expr-"=5, "legend"=0.4)

  if(sepPlots) { # matrix
     if(!missing(y))
        n <- ncol(y)
     if(!missing(segObj)) { # S4
        if(is.environment(segObj)) {
          segmentationObjectName = paste(chr, "+", sep=".")
          if(segmentationObjectName %in% ls(segObj)) {
            s <- get(segmentationObjectName, segObj)
            n <- ncol(s@y)
          }
          else { # old style list
            dat = get(paste(chr, strand, "dat", sep="."), segObj)
            n <- ncol(dat$y)
          }
        }
     }
     if(reOrder)
       ordering = seq(n,1, by=-1)
     else
       ordering = seq(1:n)
     if(n<4) {
         exprw <- exprc <- NULL
         for(i in 1:n) {
           exprw <- c(exprw, paste("expr", i, "+", sep=""))
           exprc <- c(exprc, paste("expr", i, "-", sep=""))
           }
         VPnames <- c("title", exprw, "gff+", "coord", "gff-", exprc, "legend")
         VP = c(0.8, rep(5, n), 1, 1, 1, rep(5, n), 0.4)
         names(VP) <- VPnames
      } else{
          cat("More than 4 arrays, plotting averages.\n")
          sepPlots=FALSE
      }
  }

  if(!missing(extras)) {
     indgff <-  grep("gff\\+", names(VP))
     indlegend <-  grep("legend", names(VP)) 
     VP <- c(VP[1:(indgff-1)], "extras+"=1, VP[indgff:(indlegend-1)], "extras-"=1, VP[indlegend:length(VP)])
  }
  if(!doLegend)
     VP = VP[names(VP)!="legend"]
  if(missing(gff))
     VP = VP[!(names(VP)=="gff+" | names(VP)=="gff-")]
  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
    "cp" = "#555555", "ci" = "#777777", "highlight" = "red", "threshold" = "grey")
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

  pushViewport(viewport(width=0.85, height=0.95)) ## plot margin
  pushViewport(viewport(layout=grid.layout(length(VP), 1, heights=VP)))
  for(i in 1:2) {
    strand = c("+", "-")[i]

    ## extract and treat  the data
    threshold = as.numeric(NA)

    ## Three mutually exclusive cases:
    ## 1.) y and probeAnno
    ## 2.) segObj is an environment and contains objects of S4 class "segmentation"
    ##    whose names are obtained by paste(chr, strand, sep=".")
    ## 3.) segObj is an environment and contains lists
    ##    whose names are obtained by paste(chr, strand, "dat", sep=".")
    ##
    if(!missing(y)) {
      ## case 1.
      stopifnot(is.matrix(y))
      index = get(paste(chr, strand, "index", sep="."), envir=probeAnno)
      sta   = get(paste(chr, strand, "start", sep="."), envir=probeAnno)
      end   = get(paste(chr, strand, "end",   sep="."), envir=probeAnno)
      if(!missing(extras))
        dat = list(x = (sta+end)/2,
                   y   = y[index,, drop=FALSE],
                   flag = get(paste(chr, strand, "unique", sep="."), envir=probeAnno),  
                   extras = extras[index,, drop=FALSE]) # extras not currently supported 
							# for case 2 or 3.
      else   
        dat = list(x = (sta+end)/2,
                   y = y[index,, drop=FALSE],
                   flag = get(paste(chr, strand, "unique", sep="."), envir=probeAnno))
      stopifnot(is.numeric(dat$flag))
      lengthChr = end[length(end)]
      
    } else {
      if(!is.environment(segObj))
        stop("'segObj' must be an environment.")
      
      segmentationObjectName = paste(chr, strand, sep=".")
      if(segmentationObjectName %in% ls(segObj)) {
        ## case 2: S4 class
        s = get(segmentationObjectName, segObj)
        if(!inherits(s, "segmentation"))
          stop(sprintf("'%s' must be of class'segmentation'.", segmentationObjectName))
        if(is.na(s@nrSegments))
          stop(sprintf("Slot 'nrSegments' of '%s' must not be NA.", segmentationObjectName))
        bp = s@breakpoints[[s@nrSegments]]
        dat = list(x=s@x, y=s@y, flag=s@flag, estimate = bp[, "estimate"])
        if("upper" %in% colnames(bp)) dat$upper = bp[, "upper"]
        if("lower" %in% colnames(bp)) dat$lower = bp[, "lower"]
        lengthChr <- max(s@x, na.rm=TRUE)

      } else {
        ## case 3: list 'dat' and other stuff
        dat = get(paste(chr, strand, "dat", sep="."), segObj)
        stopifnot(all(c("start", "end", "unique", "ss") %in% names(dat)))
        dat$x = (dat$start + dat$end)/2
        dat$flag = dat$unique
        lengthChr = dat$end[length(dat$end)]
        
        if("segScore" %in% ls(segObj)) {
          sgs = get("segScore", segObj)
          sgs = sgs[ sgs$chr==chr & sgs$strand==strand, c("start", "end") ]
        } else {
	  stop("This option is deprecated")
	  ##nrSegments = ...
          ##seg = get(paste(chr, strand, "seg", sep="."), segObj)
          ##th  = c(1, seg$th[nrSegments+1, 1:(nrSegments+1)])
          ##sgs = list(start  = dat$start[dat$ss][th[-length(th)]],
          ##           end    = dat$end[dat$ss][th[-1]]-1)
        }
        dat$estimate = (sgs$start[-1] + sgs$end[-length(sgs$end)]) / 2

        if("theThreshold" %in% ls(segObj))
          threshold = get("theThreshold", segObj)
      } 
    }
    ## At this point, no matter what the input to the function was,
    ##    we have the list 'dat' with elements
    ## x: x coordinate
    ## y: y coordinate
    ## flag
    ## extras (optional, for plotting an extra panel, such as p-values for each segment)
    ## estimate (optional)
    ## lower (optional)
    ## upper (optional)
    ## and possibly also a non-NA value for 'threshold'.

    ## if no region is specified, plot the whole chromosome
    if(missing(coord))
      coord = c(1, lengthChr)
    
    ## plot the data
    vpr=which(names(VP)==sprintf("expr%s", strand))
    switch(what,
      "dots" = {
      if(sepPlots) {
       ylimdata = quantile(as.vector(dat[["y"]][dat[["x"]]>=coord[1] & dat[["x"]]<=coord[2],]), c(0, 1), na.rm=TRUE)
       ylim=ylimdata 
       if(missing(ylab))
           ylab=colnames(dat$y)
       if(length(ylab)==1)
           ylab=rep(ylab, n)
       for(j in seq(1:n)) {
         datj <- dat
         datj$y <- dat$y[,ordering[j]]
         if(missing(ylab))
           ylab=colnames(dat$y)
         ## plot the data
         vpr=which(names(VP)==sprintf(paste("expr",j,"%s",sep=""), strand))
  
         plotSegmentationDots(datj, xlim=coord, ylim=ylim, ylab=ylab[ordering[j]], 
                         chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                         vpr=vpr, colors=colors, sepPlots=sepPlots,...)
         } 
      } else
         plotSegmentationDots(dat, xlim=coord, ylab=ylab, 
                         chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                         vpr=vpr, colors=colors, sepPlots=sepPlots,...)
 
       if(!missing(extras) & !missing(y)) {
             vpr2=which(names(VP)==sprintf("extras%s", strand))
             dat$y = dat$extras[,,drop=FALSE]
             plotSegmentationDots(dat, xlim=coord, chr=chr, 
                     strand=ifelse(isDirectHybe, otherStrand(strand),strand),
                     vpr=vpr2, colors=colors, colHeatmap=colExtras, 
                     ylab=ylabExtras, rowNames=rowNamesExtras,...)
        }
      },

      ## FIXME: Matt's spaghetti code needs cleanup
           
      "heatmap" = {
        plotSegmentationHeatmap(dat, xlim=coord, 
                                rowNames=rowNamesHeatmap,
                                chr=chr, 
                                strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                                vpr=vpr, colors=colors, ylab=ylab,
                                colHeatmap=colHeatmap,...)
        if(!missing(extras) & !missing(y)) {
             vpr2=which(names(VP)==sprintf("extras%s", strand))
             dat$y = dat$extras[,,drop=FALSE]
             plotSegmentationHeatmap(dat, xlim=coord, chr=chr, 
                     strand=ifelse(isDirectHybe, otherStrand(strand),strand),
                     vpr=vpr2, colors=colors, colHeatmap=colExtras, 
                     ylab=ylabExtras, rowNames=rowNamesExtras,...)
        }
      },
           stop(sprintf("Invalid value '%s' for argument 'what'", what))
    ) ## switch
  
    ## plot the features
    if(!missing(gff))
      plotFeatures(gff=gff, chr=chr, xlim=coord, strand=strand, 
                   featureExclude=featureExclude, featureColorScheme=featureColorScheme,
                   vpr=which(names(VP)==sprintf("gff%s", strand)),...) 
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
    ## this part was modified to draw arrows for transcripts rather than bars
    mt = (match(highlight$strand, c("-", "+"))-1.5)*2
    co = highlight$coord
    if(is.na(mt) || !is.numeric(co))
      stop("Invalid parameter 'highlight'.")
    strand.num <- ifelse(highlight$strand=="-",-1,1)
    grid.segments(x0=co, x1=co+(500*strand.num), y0=c(0.4,0.4)*mt, y1=c(0.4,0.4)*mt, default.units = "native", arrow=arrow(), gp=gpar(col="violetred4", lwd=4))
  }
  popViewport()

  ## title
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(VP)=="title")))
  grid.text(label=paste("Chr ", chr, sep=""), x=0.5, y=1, just="centre", gp=gpar(cex=1))
  if(!missing(main))
    grid.text(label=main, x=0.05, y=1, just="centre", gp=gpar(cex=1))
  popViewport()

  ## legend
  if(doLegend)
    plotAlongChromLegend(which(names(VP)=="legend"),
         featureColorScheme=featureColorScheme, featureExclude=featureExclude)
  
  popViewport(2)
}

## ------------------------------------------------------------
## plot Features
## ------------------------------------------------------------

plotFeatures = function(gff, chr, xlim, strand, vpr, featureColorScheme=1, featureExclude=c("chromosome", "nucleotide_match", "insertion"), featureNoLabel=c("uORF", "CDS"),...) {

  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  sel = which(gff[, "chr"] == chr &
              gff[, "strand"]  == strand &
              gff[, "start"] <= xlim[2] &
              gff[, "end"]   >= xlim[1])

  stopifnot(length(strand)==1, strand %in% c("+", "-"))
  
  ## for label, use "gene" if available, otherwise "Name"
  geneName = gff[sel, "gene"]
  featName = gff[sel, "Name"]
  featName[!is.na(geneName)] = geneName[!is.na(geneName)]

  ## split by feature type (e.g. CDS, ncRNA)
  feature  = as.character(gff[sel, "feature"])
  featsp = split(seq(along=sel), feature)

  ## There are now five different cases, and we need to deal with them:
  ## - ignorable features, given by featureExclude
  ## - genes: a horizontal line + name
  ## - introns: a caret
  ## - CDS: a box + no name
  ## - all others: a colored box + name

  ## in this vector we save those features for which we want to have names
  whnames = integer(0)

  ## 1. drop the ignorable ones
  featsp = featsp[ ! (names(featsp) %in% featureExclude) ]
  
  ## 2. gene: just a horizontal line + name
  wh = ("gene" == names(featsp))
  if(any(wh)) {
    i = featsp[["gene"]]
    s = sel[i]
    grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, y1 = 0,
                  default.units = "native", gp = gpar(col="#a0a0a0"))
    whnames = i
    featsp = featsp[!wh]
  }

  ## 3.introns
  wh = ("intron" == names(featsp))
  if(any(wh)) {
    i = featsp[["intron"]]
    s = sel[i]
    mid = (gff$start[s]+gff$end[s])/2
    wid = (gff$end[s]-gff$start[s])/2 
    for(z in c(-1,1))
      grid.segments(x0 = mid,
                    x1 = mid+z*wid,
                    y0 = 1.20*c("+"=1, "-"=-1)[strand],  ## istrand is 1 or 2
                    y1 = 0.95*c("+"=1, "-"=-1)[strand],
                    default.units = "native",
                    gp = gpar(col="black"))
     featsp = featsp[!wh]
  } ## if
  
  ## 4. colors for boxes
  ## check that we know how deal with all features
  featCols = featureColors(featureColorScheme)

  whm = names(featsp) %in% rownames(featCols)
  if(!all(whm))
    warning("Don't know how to handle feature of type(s) '", paste(names(featsp)[!whm], collapse=", "), "' in gff.", sep="")

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
  darken = function(x, factor=0.5) {
    wh = which(!is.na(x))

    hex = sapply(x[wh], substring, first=c(2,4,6), last=c(3,5,7))
    hex = apply(hex, 2, function(h) as.integer(factor*as.integer(paste("0x", h, sep=""))))

    res = rep(as.character(NA), length(x))
    res[wh] = apply(hex, 2, function(h) sprintf("#%02x%02x%02x", h[1], h[2], h[3]))
    return(res)
  }
  
  border = ifelse(darkenborder, darken(fill), fill)
  
  res = data.frame(fill=I(fill),
    col =I(border))
  rownames(res)=names(defaultColors) 
  return(res)
}
