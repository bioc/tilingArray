library(tilingArray)
options(error=recover)

sqr = function(x){x*x}


## what would it cost to move the changepoints?
wiggle = function(y, cp, J, alpha=0.95, double.check=TRUE) {
  
  ## create of a list whose elements are the data values within each segment
  cp1 = c(1, cp)
  segY = lapply(seq(along=cp), function(k) {
    y[cp1[k]:(cp1[k+1]-1)]
  })
  segSum = sapply(segY, sum)
  segSS  = mapply(function(x,s) sum(sqr(x-s/length(x))), segY, segSum)
  segLen = listLen(segY)
  
  ## J = -log( SS / n) !
  stopifnot(abs(J + log(sum(segSS)/length(y))) < 1e-10)
  
  ## estimated variance (unbiased version)
  s2 = sum(segSS) / (length(y)-length(cp))
  
  ## By how much can the sum of squares increase?
  ## Denote by L the negative log-likelihood, by SS the sum of
  ## squared residuals, and by n the number of residuals, then
  ##   Delta L  = n/2 log (Delta SS/n)   <=>
  ##   Delta SS = n * exp(2L/n)
  ##
  ## We assume that the likelihood ratio statistic follows
  ##   -2 log(L_rest/L_abs) ~ Chi(1)^2
  ## implying a critical value for -log(L_rest/L_abs)
  ##   of qchisq(alpha, df=1)/2
  ## Hence the critical value for SS is
  ##   n * exp(qchisq(alpha, df=1)/n)
  
  ## qchisq(alpha, df=1)
 
  dsl = dsr = vector(mode="list", length=length(cp)-1)
 
  for(j in 1:length(dsl)) {

    th  = cp[j]
    thl = if(j>=2) {cp[j-1]} else {1}
    thr = cp[j+1]
    
    yj  = y[thl:(thr-1)]
    sj  = (segSS[j]+segSS[j+1])
    
    rl = numeric(th-thl)
    for(i in 1:(th-thl)) {
      ## extend to the left by i
      yl  = sum(y[th-(1:i)])
      m1l = (segSum[j]   - yl) / (segLen[j]  -i)
      m2l = (segSum[j+1] + yl) / (segLen[j+1]+i)
      nml = rep(c(m1l, m2l), c(th-thl-i, thr-th+i))
      ssl = sum(sqr( yj-nml ))
      rl[i] = ssl - sj
    }
    stopifnot(all(rl>=0))
    dsl[[j]] = rl
    
    rr = numeric(thr-th)
    for(i in 1:(thr-th)) {
      ## extend to the right by i
      yr  = sum(y[th+(0:(i-1))])
      m1r = (segSum[j]   + yr) / (segLen[j]   +i)
      m2r = (segSum[j+1] - yr) / (segLen[j+1] -i)
      nmr = rep(c(m1r, m2r), c(th-thl+i, thr-th-i))
      ssr = sum(sqr( yj-nmr ))
      rr[i] = ssr - sj
    }
    stopifnot(all(rr>=0))
    dsr[[j]] = rr

    if(double.check) {
      ## double-check (slow - can be removed eventually)
      for(i in c(-seq(along=logical(length(rl)-1)),
                 seq(along=logical(length(rr)-1)))) {
        cp2 = cp1
        cp2[j+1] = cp2[j+1]+i
        segY2 = lapply(seq(along=cp), function(k) {
          y[cp2[k]:(cp2[k+1]-1)]
        })
        segSum2 = sapply(segY2, sum)
        segSS2  = mapply(function(x,s) sum(sqr(x-s/length(x))), segY2, segSum2)
        stopifnot(abs(sum(segSS2) - sum(segSS) - ifelse(i>0, rr[i], rl[-i])) < 1e-10)
      } ## for i
    } ## if double.check
  } ## for j
  return(list(dsl=dsl, dsr=dsr))
}

simulate = function(n=100, wid=20, sigma=1) {
  repeat {
    y = rnorm(n, mean=0, sd=sigma)
    fs = findSegments(y, maxcp=2, maxk=n)
    estcp = fs$th[2,]
    if((estcp[1]>wid) && (estcp[2]-estcp[1]>=wid))
      break
  } 
  z = wiggle(y, cp=estcp, J=fs$J[2], double.check=!TRUE)
}


nrep = 1000
n = 200
kost=matrix(as.numeric(NA), nrow=2*n+1, ncol=nrep)
for(i in 1:nrep) {
  cat(i, "")
  s = simulate(n=n)
  dsl = s$dsl[[1]]
  dsr = s$dsr[[1]]
  kost[n:(n-length(dsl)+1), i] = dsl
  kost[n+1+(1:length(dsr)), i] = dsr
}

par(mfrow=c(1,2))
plot((rowMeans(kost,na.rm=T)))
image(kost)

stop()

par(mfrow=c(2,1))
plot(y)
abline(v=estcp, col="red")

mapply(function(dsl, dsr){
  plot(seq(-length(dsl), length(dsr)),
       c(rev(dsl), 0, dsr), type="b", pch=16, col="blue", ylim=c(0,20))
}, z$dsl, z$dsr)


##l = c(10,  17, 3, 27, 10, 13)
##m = c( 2, 2.2, 4, 0,   2, 1.4)
##S = length(l)
##maxk = 30

##cp = 1+cumsum(l)
##s  = 0.6

##y = unlist(lapply(seq(along=l), function(i) {
##  rnorm(l[i], mean=m[i], sd=s)
##}))

