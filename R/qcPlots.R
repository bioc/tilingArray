myWrite = function(x, nm, imgdir=NULL, pgm=FALSE) {
  if(pgm) {
    f1 = file.path(ifelse(is.null(imgdir), ".", imgdir), paste(nm, ".pgm", sep=""))
    write.pnm(x, file=f1, maxval=255)
  }
  else { ## convert from pgm to .jpg
    f1 = paste(tempfile(), ".pgm", sep="")
    f2 = file.path(ifelse(is.null(imgdir), ".", imgdir), paste(nm, ".jpg", sep=""))
    write.pnm(x, file=f1, maxval=255)
    cmd = paste("convert", f1, "-quality 100", f2, "&")
    system(cmd)
  }
}

myPixmap = function(z, rz=range(z), nr = 2560, nc = 2560) {
  nz  = (z-rz[1])/diff(rz)
  nz[is.na(nz)] = 0
  pixmapGrey(data=matrix(nz, nrow=nr, ncol=nc, byrow=TRUE))
}

matDensities <- function(X) {
        densXY <- function(Z) {
            zd <- density(Z, na.rm = TRUE)
            x <- zd$x
            y <- zd$y
            cbind(x, y)
        }
        out <- apply(X, 2, densXY)
        outx <- out[(1:(nrow(out)/2)), ]
        outy <- out[(((nrow(out)/2) + 1):nrow(out)), ]
        list(X = outx, Y = outy)
    }

qcPlots <- function(x, html=TRUE, plotdir=NULL, probeAnno, gff, 
                    chr=4, coord=c(230000,245000),
  		    nr = 2560, nc = 2560,
		    ylimchrom=c(5,16), nucleicAcid, pmindex, pgm=TRUE,
                    ext=".cel", ranks=FALSE, ...) {

    cat("Generating summary plots\n")
    narrays = nrow(pData(x))
    l2x  = log2(exprs(x))

    # boxplots of log2 intensities for PM probes
    cat("\tBoxplots of PM intensities")
    if(missing(nucleicAcid)) {
        nucleicAcid = NULL
      if(!is.null(x$nucleicAcid))
        nucleicAcid = x$nucleicAcid
      if(!is.null(x$NucleicAcid))
        nucleicAcid = x$NucleicAcid
      if(is.null(nucleicAcid))
        nucleicAcid = rep("RNA", narrays)    
      }
    if(missing(pmindex))
      pmindex = PMindex(probeAnno)
    l2xPM = l2x[pmindex,, drop=FALSE]
    z   = lapply(1:ncol(l2xPM), function(i) l2xPM[,i])
    colors = brewer.pal(9, "Set1")[as.integer(factor(nucleicAcid))]
    filename = file.path(ifelse(is.null(plotdir), ".", plotdir), "boxplot.png")
    png(filename, width=640, height=480)
    boxplot(z, col=colors, outline=FALSE, las=2, main = "Boxplots of PM intensities", ylab=expression(log[2](Intensity)))
    dev.off()
    rm(z)
    cat(" .... complete\n")
 
    # density plots of log2 intensities for PM probes
    cat("\tDensity plot of PM intensities")
    dens.x <- matDensities(l2xPM)
    filename = file.path(ifelse(is.null(plotdir), ".", plotdir), "densities.png")
    png(filename, width=640, height=480)
    matplot(dens.x$X, dens.x$Y, xlab=expression(log[2](Intensity)),
            ylab = "Density", col=colors,
            main = "Histogram of PM intensities",
            type = "l", lwd = 2, lty = 1)
    dev.off()
    rm(l2xPM)
    cat(" .... complete\n")

    # individual image, density and along chromosome dot plots (for region specified by user)
    cat("\nGenerating per array plots\n")
    rlx = range(l2x)
    lx  = (l2x-rlx[1])/diff(rlx)

    files <- NULL
    if(!is.null(x$file))
      files = x$file
    if(!is.null(x$File))
      files = x$File
    if(is.null(files))
      files = paste("array", seq(1:narrays), ext, sep="")

    for(i in 1:narrays) {
      cat(files[i], ":\nimageplot\t")
      filename = gsub(" ", "_", files[i])
      if(ranks)
      myWrite(myPixmap(rank(exprs(x)[,i])), sub(ext, ".rank", filename),
              imgdir=plotdir, pgm=pgm)
      if(!ranks)
      myWrite(myPixmap(lx[,i]),  sub(ext, ".log", filename),
              imgdir=plotdir, pgm=pgm)
      cat("density plot\t")     
      filename = file.path(ifelse(is.null(plotdir), ".", plotdir), 
                  paste(sub(ext, ".density", files[i]), ".png", sep=""))
      png(filename)
      if(narrays>1)
         matplot(dens.x$X[,i], dens.x$Y[,i], xlab=expression(log[2](Intensity)),            
                  ylab = "Density", col=colors[i], main = "Histogram of PM intensities", 
                  type = "l", lwd = 2, lty = 1)
      else
         matplot(dens.x$X, dens.x$Y, xlab=expression(log[2](Intensity)), ylab = "Density",
                  col=colors[i], main = "Histogram of PM intensities", type = "l", 
                  lwd = 2, lty = 1)
      dev.off()
      cat("along chromosome plot")     
      filename = file.path(ifelse(is.null(plotdir), ".", plotdir), paste(sub(ext, ".gencoord", files[i]), ".jpg", sep=""))
      jpeg(filename, width=960, height=480)
      plotAlongChrom(what="dots", chr=chr, coord=coord, y=matrix(l2x[,i], nc*nr,1),
           probeAnno=probeAnno, gff=gff, ylim=ylimchrom, ...)
      dev.off()
      cat(" .... complete\n")
      }
   rm(l2x)
   rm(lx)
   rm(rlx)
   gc()

   if(html) {
    cat("Creating html summary page")
     html <- list(NULL)

     html[[1]] = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n<title>Quality assessment of tiling array data</title>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n</head>"

     html[[2]] = "<h1>Quality assessment of tiling array data</h1>\n"
     html[[3]] = paste("<p>",date(),"</p>\n", sep="")
     html[[4]] = "<h2>Summary plots</h2>\n"
     html[[5]] = "<img src=\"boxplot.png\"><br>\n" 
     html[[6]] = "<img src=\"densities.png\"><br>\n"
     html[[7]] = "<h2>Individual plots</h2>\n"
     if(!is.null(x$SampleID) & !is.null(x$Strain)) 
          html[[8]] = "<table border=1>\n<tr bgcolor=\"lightgray\"><td><center><strong>Array</strong></center></td><td><center><strong>Filename</strong></center></td><td><strong>Strain</strong></td><td><strong>SampleID</strong></td><td><strong>Plots</strong></td></tr>\n"
     else
        html[[8]] = "<table border=1>\n<tr bgcolor=\"lightgray\"><td><center><strong>Array</strong></center></td><td><center><strong>Filename</strong></center></td><td><center><strong>Plots</strong></center></td></tr>\n"
     for(i in 1:narrays) {
        if(!is.null(x$SampleID) & !is.null(x$Strain)) 
          html[[i+8]] = paste("<tr><td><center>", i, "</center></td><td>", files[i], "</td><td>", x$Strain[i], "</td><td>", x$SampleID[i], "</td><td>
        <a href=\"", sub(ext, ".log", files[i]), ".jpg\">image</a>,", "<a href=\"", sub(ext, ".density", files[i]), ".png\">", "density</a>,", "<a href=\"", sub      (ext, ".gencoord", files[i]), ".jpg\">", "alongchrom</a></td></tr>\n", sep="")
          
        else 
          html[[i+8]] = paste("<tr><td><center>", i, "</center></td><td>", files[i], "</td><td> <a href=\"", sub(ext, ".log", files[i]), ".jpg\">image</a>,", "<a href=\"", sub(ext, ".density", files[i]), ".png\">", "density</a>,", "<a href=\"", sub      (ext, ".gencoord", files[i]), ".jpg\">", "alongchrom</a></td></tr>\n", sep="")
        }
     html[[i+9]] = "</table>\n</body>\n</html>"

     if(pgm)
       html = gsub("log.jpg", "log.pgm", html)
     if(ranks)
       html = gsub("log.", "rank.", html)
     filename = file.path(ifelse(is.null(plotdir), ".", plotdir), "qcsummary.htm")
     writeLines(unlist(html), con=filename)    
     cat(" .... complete\nOpen qcsummary.htm to view results.\n")
    }
}
