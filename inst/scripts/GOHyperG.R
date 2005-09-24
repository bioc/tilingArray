GOHyperG = function(candidates) {
  require("GO")
  require("annotate")
  
  if(!("Ontology_term" %in% names(gff)))
    gff$Ontology_term=getAttributeField(gff[, "attributes"], "Ontology_term")
  
  ctrls  = gff[ gff[, "feature"]=="gene", "Name"]
  ctrls  = unlist(strsplit(ctrls, split=", "))
  
  
  e = new.env(hash=TRUE)
  for(j in ls(GOMFANCESTOR))
    assign(j, get(j, GOMFANCESTOR), envir=e)
  for(j in ls(GOBPANCESTOR))
    assign(j, get(j, GOBPANCESTOR), envir=e)
  for(j in ls(GOCCANCESTOR))
    assign(j, get(j, GOCCANCESTOR), envir=e)
  stopifnot(length(ls(e))==length(ls(GOMFANCESTOR))+
            length(ls(GOBPANCESTOR))+length(ls(GOCCANCESTOR)))
  
  ## for each gene in 'x', get the GO classes
  whg = which(gff[, "feature"]=="gene")
  getGO = function(x) {
    mt  = match(x, gff[whg, "Name"])
    rv = strsplit(gff[whg[mt], "Ontology_term"], split=",")
    stopifnot(!any(sapply(rv, function(x) any(duplicated(x)))))
    
    ## extend by ancestors
    rv = sapply(rv, function(v) {
      if(any(is.na(v))) {
        stopifnot(length(v)==1)
        k = character(0)
      } else {
        k  = mget(v, e, ifnotfound=list(character(0)))
        k  = sort(unique(unlist(k)))
      }
      return(k)
    })
    rv
  }
  
  ## if(!exists("ctGO")){
    asGO = getGO(candidates)
    ctGO = getGO(ctrls)
  ##}

  stopifnot(!any(sapply(asGO, function(z) any(duplicated(z)))), 
            !any(sapply(ctGO, function(z) any(duplicated(z)))))
  
  asTab = table(unlist(asGO))
  ctTab = table(unlist(ctGO))
  
  allGO = sort(unique(c(unlist(asGO), unlist(ctGO))))
  nr = matrix(NA, nrow=length(allGO), ncol=2)
  rownames(nr) = allGO
  
  nr[names(asTab), 1] = asTab
  nr[names(ctTab), 2] = ctTab
  
  ## white balls: genes with GO Term
  ## black balls: genes without GO Term
  ## phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE):
  ##     m: the number of white balls in the urn.  : nr[,2]
  ##     n: the number of black balls in the urn.  : length(ctGO)-nr[,2]
  ##     k: the number of balls drawn from the urn.: length(asGO)
  ph = phyper(nr[,1], m=nr[,2], n=length(ctGO)-nr[,2], k=length(asGO), lower.tail=FALSE)
  
  thSlop = length(candidates) / length(ctrls)
  
  ksel = which(nr[,1]>=4 & nr[,1]<length(asGO))
  ksel = ksel[order(ph[ksel])[1:40]]
  
  for(k in ksel) {
    L = rownames(nr)[k]
    cat("------------------------------------------------------------\n")
    print(get(L, GOTERM))
    cat("\nIn genome: ", nr[L,2], ", expected: ", round(nr[L,2]*thSlop,1), ", found: ", nr[L, 1],
        ", p=", format.pval(ph[L]), "\n",
        "Genes: ",
        paste(sort(replaceSystematicByCommonName(candidates[sapply(asGO, function(x) L %in% x)])), collapse=" "),
        "\n\n", sep="")
    
  }
  
  par(mfrow=c(1,2))
  cols=rep("grey", nrow(nr))
  cols[ksel]="orange"

  xmax = c(max(nr[,1], na.rm=TRUE),  20)
  for (i in seq(along=xmax)) {
    plot(nr, pch=16, xlab="antisense genes", ylab="control: all genes",
         main=c("Frequency of occurence of GO classes", "zoom in")[i], 
         xlim=c(0, xmax[i]), ylim=c(0, xmax/thSlop), col=cols)
    abline(a=0, b=1/thSlop, col="blue")
  }
}
