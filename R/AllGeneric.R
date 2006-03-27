##
## S3 generic functions "confint" and "logLik" are defined in package "stats",
##   here we try to keep consistent with that:
##
## > confint
## function (object, parm, level = 0.95, ...)
## UseMethod("confint")
## <environment: namespace:stats>
##
## > logLik
## function (object, ...)
## UseMethod("logLik") 
## <environment: namespace:stats>
##

if(!isGeneric("confint"))
  setGeneric("confint", useAsDefault=confint)
##  setGeneric("confint", function(object, parm, level = 0.95, ...) standardGeneric("confint"))

if(!isGeneric("logLik"))
  setGeneric("logLik", useAsDefault=logLik)
##  setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

if(!isGeneric("plot"))
  setGeneric("plot", useAsDefault=plot)

if(!isGeneric("show"))
  setGeneric("show", useAsDefault=show)
