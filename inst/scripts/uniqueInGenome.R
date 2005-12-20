if(!exists("probeAnno"))
  load("probeAnno.rda")
library("Biobase")

chrstr = paste(rep(1:16, each=2), rep(c("+","-"), 16), sep=".")
uni = mget(paste(chrstr, "unique", sep="."), probeAnno)
sta = mget(paste(chrstr, "start", sep="."), probeAnno)
end = mget(paste(chrstr, "end", sep="."), probeAnno)
stopifnot(all(listLen(uni)==listLen(sta)), all(listLen(uni)==listLen(end)))

tot = 0
dup = 0
for(i in seq(along=uni)) {
  y = (uni[[i]]>1)
  x = (sta[[i]]+end[[i]])/2
  app = approx(x, y, xout=seq(min(x), max(x), by=1), rule=2)$y
  s   = sum(app)
  tot = tot + length(app)
  dup = dup + s
  cat(sprintf("%s : %d of %d = %3.1f percent\n", chrstr[i], as.integer(s), length(app), s/length(app)*100))
}

cat(sprintf("total : %d of %d = %3.1f percent\n", as.integer(dup), as.integer(tot), dup/tot*100))
