options(error=recover, warn=2)
library("tilingArray")
library("prada")

indir = "segmentation-050209v4"
js = paste(rep(1:16, each=2), rep(c("+", "-"), 16), "rda", sep=".")

if(!exists("s")) {
  s  = new.env()
  cat("Loading.\n")
  for(j in js) {
    load(file.path(indir, j))
    assign(sub("rda", "seg", j), seg, envir=s)
    assign(sub("rda", "dat", j), dat, envir=s)
  }
}


nrBasePerSeg = 1500
j = js[7]

dat = get(sub("rda", "dat", j), envir=s)
seg = get(sub("rda", "seg", j), envir=s)

cp = round(max(dat$x)/nrBasePerSeg)
th = c(1, seg$th[cp, 1:cp])
ct = cut(1:length(dat$y), breaks=th-1)
sp = split(dat$y, ct)
v  = sapply(sp, var)
m  = sapply(sp, mean)
n  = listLen(sp)

## t-statistic:
## (m1-m2) / sqrt((var1+var2)*(n1+n2)/((n1+n2-2)*n1*n2))

t = numeric(length(sp))
for (z in 2:length(sp)) {
  s = sqrt(  ((n[z-1]-1)*v[z-1]+(n[z]-1)*v[z]) / (n[z]+n[z-1]-2) )
  t[z] = (m[z]-m[z-1]) / (s*sqrt(1/n[z]+1/n[z-1]))
}

z = sample(2:length(sp), 1)
stopifnot(abs(t[z] - t.test(sp[[z]], sp[[z-1]], var.equal=TRUE)$statistic) < 1e-9)



if(FALSE) {
  v = sapply(js, function(j)
     sd(get(sub("rda", "dat", j), envir=s)$y))
}

## see whether variance depends on position along chromosome
## (it doesn't seem to)
if(FALSE) {
  par(mfrow=c(4,4), mai=c(0.3, 0.4, 0.4, 0.01))
  for(j in js[17:32]) {
    dat = get(sub("rda", "dat", j), envir=s)
    dy  = diff(dat$y)
    dy  = dy[abs(dy)<=2]
    smoothScatter(seq(along=dy), dy, pch=".", main=j)
  }
}
