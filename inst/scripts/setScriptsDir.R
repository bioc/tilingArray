scriptsDir = function(x) {
  whoami = system("whoami", intern=TRUE)
  prefix = c("huber"="scripts", "toedling"="jscripts")[whoami]
  if(is.na(prefix))
    stop(paste("Don't know how to deal with '", whoami, "'", sep=""))
  file.path(prefix, x)
}

functionsDir = function(x) {
  whoami = system("whoami", intern=TRUE)
  prefix = c("huber"    = "/ebi/research/huber/users/huber/madman/Rpacks/tilingArray/R",
             "toedling" = "/ebi/research/huber/users/joern/tilingArray/R")[whoami]
  if(is.na(prefix))
    stop(paste("Don't know how to deal with '", whoami, "'", sep=""))
  file.path(prefix, x)
}
