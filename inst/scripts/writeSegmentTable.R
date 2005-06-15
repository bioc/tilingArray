
alongChromWidth = 15e3
alongChromStep  =  5e3

## this function maps a chromosome number and start and end coordinates
## to a file name
mapCoord2Plot = function (chr, start, end) {
  mid  = (start+end-alongChromWidth)/2
  mid[mid<0]=0
  pst  = as.integer(alongChromStep/1e3*ceiling(mid/alongChromStep))
  sprintf("%02d_%04d", as.integer(chr), pst)
}

replaceSystematicByCommonName = function(z) {
  sapply(strsplit(z, split=", "), function(x) {
    mt = match(x, gff[, "Name"])
    gene = gff[mt, "gene"]
    paste(ifelse(is.na(gene)|gene=="", x, gene), collapse=", ")
  })
}

writeSegmentTable = function(sgs, title, fn, sortBy, sortDecreasing=FALSE,
  colors=c("#d0d0d0", "#ffffff"), interact, HTML=TRUE) {

  fn = paste(fn, ifelse(HTML, "html", "txt"), sep=".")
  
  if(nrow(sgs)==0)
    stop("'sgs' has 0 rows")
  if(interact)
    cat("Writing ", fn, " (", nrow(sgs), " rows)\n", sep="")
  
  for(cn in c("featureInSegment", "mostOfFeatureInSegment", "overlappingFeature", "oppositeFeature"))
    sgs[, cn] =  replaceSystematicByCommonName(sgs[, cn])
  
  rowOrd = switch(sortBy,
    "category-level" = {
      stopifnot(all(c("category","level") %in% names(sgs)))
      order(sgs[, "category"], -sgs[, "level"], decreasing=FALSE)
    },
    "goodUTR" = {
      stopifnot("goodUTR" %in% names(sgs))
      order(sgs[, "goodUTR"])
    },
    stop("Zapperlot"))
    
  colOrd = c("segID", "category", "overlap", "overlappingFeature", "oppositeFeature",
    "chr", "strand", "start", "end", "length", "level", "utr5", "utr3",
    "zLeft", "zRight",
    "featureInSegment", "mostOfFeatureInSegment", "overlapFeatAll", 
    "oppositeExpression", "distLeft", "distRight", "simpleCatg", 
    "frac.dup")
  colOrd = c(sortBy, colOrd[!(colOrd %in% sortBy)])

  inOrd = colnames(sgs) %in% colOrd
  inSgs = colOrd %in% colnames(sgs)

  if(interact){
    if(!all(inOrd))
      cat("Warning: column(s)", colnames(sgs)[!inOrd], "not in 'colOrd'.\n")
    if(!all(inSgs))
      cat("Warning: column(s)", colOrd[!inSgs], "not in the table.\n")
  }
  sgs = sgs[rowOrd, colOrd[inSgs]]
  
  out = matrix(as.character(NA), nrow=nrow(sgs), ncol=1+ncol(sgs))
  colnames(out) = c("plot", colnames(sgs))
  
  for(i in 1:ncol(sgs)) {
    if(is.numeric(sgs[[i]])) {
      z = sgs[[i]]
      if(any(abs(round(z)-z)>0.01, na.rm=TRUE))
        out[, i+1] = paste(signif(z, 3))
      else
        out[, i+1] = paste(as.integer(z))        
    } else {
      out[, i+1] = as.character(sgs[[i]])
    }
  }

  ## create hyperlink
  hyper = mapCoord2Plot(sgs$chr, sgs$start, sgs$end)
  hyper = paste("<a href=\"", hyper, ".jpg\" target=\"", hyper, "\">jpg</a> ", 
    "<a href=\"", hyper, ".pdf\" target=\"", hyper, "\">pdf</a>", sep="")
  out[, 1] = hyper
    
  out = t(out)
  
  con=file(fn, open="wt")
  if(HTML) {
    cat("<html><STYLE>", 
        "<!--TD { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 14px;}-->",
        "<!--H1 { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 22px;}-->",
        "</STYLE>", "<head>", paste("<TITLE>", title, "</TITLE>", sep=""),
        "</head>", "<body bgcolor=#ffffff>", file = con, sep = "\n")
          if (title!="") 
            cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n", 
                file = con, sep = "\n")
    cat("<CENTER>\n<TABLE BORDER=0>\n", file = con)
    cat("<TR>\n", file=con)
    cat(paste("<TD><B><i>", rownames(out) , "</i></B></TD>", collapse=""), sep="", file=con)
    cat("</TR>\n", file=con)
    for(i in 1:ncol(out)) {
      cat("<TR BGCOLOR=\"", colors[(i-1) %% length(colors) + 1], "\">\n", sep="", file=con)
      cat(paste("<TD>", out[,i] , "</TD>", collapse="", sep=""), file=con)
      cat("</TR>\n", file=con)
    }
    
    cat("</CENTER>\n</TABLE>\n", file = con)
    cat("<hr><i>Wolfgang Huber</i> -- ", date(), "</html>", file = con)
  } else {
    cat(paste(rownames(out) , collapse="\t"), "\n", sep="", file = con)
    for(i in 1:ncol(out))
      cat(paste(out[,i] , collapse="\t"), "\n", sep="", file = con)
  }
  close(con)
}
