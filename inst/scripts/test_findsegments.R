library("tilingArray")
options(error=recover)
options(warn=0)



x = NULL
nseg = round(runif(1)*4+2)
b = 0
for(i in 1:nseg) {
  b      = b + (runif(1)+.5)*sign(runif(1)-0.5)
  lenseg = round(runif(1)*7+3)
  x      = c(x, b + rnorm(lenseg, sd=0.1))
}
plot(x)



maxk = 11
Km   = 7

seg  = findsegments(x, maxk=maxk, Km=Km)
abline(v=seg$th[nseg,], col="red")

## plot(Y)




###################################################################
## IN THE FOLLOWING SOME FUNCTIONS THAT IMPLEMENT THE STUFF IN R ##
###################################################################
 
if(FALSE){
##------------------------------------------------------------
## which.min with better handling of vectors with all Inf
##------------------------------------------------------------
my.which.min = function(...) {
  wh = which.min(...)
  if(length(wh)==0)
    wh=1
  return(wh)
}

##------------------------------------------------------------
## th  : coordonnees de la fin des ruptures
## G   : matrice de coût
## Km  : nb de segments max que l'on veut
##------------------------------------------------------------
dynprog = function (fG, n, Km){
  mI =  matrix(Inf, nrow=Km,   ncol=n)	
  mt  = matrix(0,   nrow=Km-1, ncol=n) 
  mI[1,] = fG(1,1:n)
  
  for (k in 2:(Km-1)) {	## k - the number of segments		
    for (L in k:n) {    ## L  - candidate jump point
      z         = mI[k-1, 1:(L-1)] + fG(2:L, L)
      mI[k,L]   = min(z)
      mt[k-1,L] = my.which.min(z)
    }
  }
  z          = mI[Km-1, 1:(n-1)] + fG(2:n,n)
  mI[Km,n]   = min(z)
  mt[Km-1,n] = my.which.min(z)
  J = mI[,n]

  ## for (i in 1:n) print(fG(i, 1:n))
  ## print(signif(mI, 3))

  ## matplot(t(mI[1:5,])) ## looks good
  
  ## Compute the change-points instants ***
  th = diag(rep(n,Km))
  for (K in 2:Km) {
    for (k in seq(K-1, 1, by=-1)) {
      th[K,k] = mt[k,th[K,k+1]]
    }
  }
  return(list(J=J,th=th))
}

## *****************************************************************************
## This function realizes the segmentation according to the model
## G is the matrix containing the likelihood for every possible segment (cost matrix)
## J is the contrast (related to the likelihood) for K=1...Km
## th contains the breakpoint coordinates for K=1...Km
## lmin is the minimum length of a segment (lmin>1 if model=1)
## *****************************************************************************
test_findsegments = function(Y, Km, maxk) {
  G = Gmean(Y, maxk=maxk)

  fG = function(i, j) {
    stopifnot(all(i>=1&i<=ncol(G)))
    k  = j-i+1
    ik = nrow(G)*(i-1) + k
    n  = length(k)
    rv = rep(Inf, n)
    flag = (k <= nrow(G))
    rv[flag] = G[ik[flag]]
    return(rv)
  }
  
  n  = length(Y)
  rv = dynprog(fG, n, Km) ## list with elements "J" and "th"
  rv$J=log(rv$J/n)
  return(rv)
}


}
