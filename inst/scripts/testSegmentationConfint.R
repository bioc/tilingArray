library(tilingArray)
options(error=recover)

for(f in dir("~/madman/Rpacks/tilingArray/R", pattern=".R$", full.names=TRUE))source(f)

set.seed(123)

lseg  = 12
steps = seq(1, 4, length=4)
steps = steps * sign(runif(length(steps))-0.5)
dat   = rep(cumsum(steps), each=lseg) + rnorm(lseg*length(steps))
dat   = cbind(dat)
##dat <- matrix(c(rnorm(10,0,1), rnorm(20,2,1), rnorm(5,0.5,0.5),
##                rnorm(10,1,1), rnorm(20,2,1)), ncol=1)

par(mfrow=c(2,1))

seg1 = segment(dat, maxseg=10, maxk=nrow(dat))
seg2 = confint(seg1, 2:5)

plot(seg2, 5)

#debug(confint.segmentation)
old1 = findSegments(dat, maxcp=10, maxk=nrow(dat))
old2 = confint.segmentation(old1, 2:5)
plot.segmentation(old2, 5)

