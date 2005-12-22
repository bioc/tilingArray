##
## An S3 generic function "confint" is defined in package "stats",
##   here we try to keep consistent with that:
##
## function (object, parm, level = 0.95, ...)
## UseMethod("confint")
## <environment: namespace:stats>
##
##
setGeneric("confint", function(object, parm, level = 0.95, ...) standardGeneric("confint"))


