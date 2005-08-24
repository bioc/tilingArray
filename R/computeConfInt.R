# function to compute confidence intervals for Picard's segmentation
#  it's a modified version of 'confint.breakpointsfull' from package strucchange

computeConfInt <- function (object, parm = NULL, level = 0.95, breaks = NULL, het.reg = FALSE, het.err = TRUE, ...)
{
  # check arguments:
  stopifnot(all(c("y","breakpoints","residuals") %in% names(object)))
  
  #X <- object$X  # vector of 1s, maybe offset in formula
  y <- as.matrix(object$y)  # numeric values to do segementation on, in our case:
                 #  intensities along genome
  n <- nrow(y) #object$nobs # length of values to segment (no. of intensities)
  X <- matrix(1, nrow=n,ncol=1)
  a2 <- (1 - level)/2 
  if (!is.null(parm) & !is.null(breaks))
    warning("`parm' and `breaks' are both specified: `breaks' is used")
  else if (!is.null(parm))
    breaks <- parm
  myfun <- function(x, level = 0.975, xi = 1, phi1 = 1, phi2 = 1)
    (pargmaxV(x, xi = xi, phi1 = phi1, phi2 = phi2) - level)
  myprod <- function(delta, mat)
    as.vector(crossprod(delta,mat) %*% delta)

  bp    <- object$breakpoints
  nbp   <- length(bp) # number of breakpoints
  upper <- rep(0, nbp)
  lower <- rep(0, nbp)
  bp    <- c(0, bp, (n+1)) # 0 and (n+1) are also considered breakpoints
                           #  changed also because of index denoting start of
                           #  next seg, not end of this one.

  res   <- as.matrix(object$residuals)
  #sigma1 <- sigma2 <- sum(res^2)/n  # sum divided by n
  sigma1 <- sigma2 <- mean(crossprod(res^2, X)/n)  # sum divided by n

  #Q1 <- Q2 <- crossprod(X)/n # amounts to 1 in the example
  Q1 <- Q2 <- 1
 
  # removed old construction with 'vcov'
  Omega1 <- Omega2 <- sigma1 * Q1

  xi <- 1
  X2 <- X[(bp[1] + 1):bp[2], , drop = FALSE]
  y2 <- y[(bp[1] + 1):bp[2]]
  fm2 <- lm(y2 ~ 0 + X2)
  beta2 <- coef(fm2) # mean of first segment

  if (het.reg)  # heterogeneous distribution of segment regressors?
    Q2 <- crossprod(X2)/nrow(X2)
  if (het.err){   # heterogeneous distribution of segment errors?
    sigma2 <- sum(residuals(fm2)^2)/nrow(X2)
    Omega2 <- sigma2 * Q2
  }
 
  for (i in 2:(nbp + 1)) { # for all segment breakpoints

    X1 <- X2
    y1 <- y2
    beta1 <- beta2
    sigma1 <- sigma2
    Q1 <- Q2
    Omega1 <- Omega2

    # get values of the current segment:
    X2 <- X[bp[i]:(bp[i + 1]-1), , drop = FALSE]
    #X2 <- 1
   
    y2 <- y[bp[i]:(bp[i + 1]-1)]
    # old version changed because our index denotes already start of next seg

    fm2 <- lm(y2 ~ 0 + X2) # just computes the mean
    beta2 <- coef(fm2)     #  this is the mean
    delta <- beta2 - beta1 #  this is the difference of means between the segs
    if (het.reg)
      Q2 <- crossprod(X2)/nrow(X2)    # this is 1
    if (het.err) {
      sigma2 <- sum(residuals(fm2)^2)/nrow(X2) # average square residual in seg
      Omega2 <- sigma2 * Q2   # equal to sigma2
    }
    Oprod1 <- myprod(delta, Omega1)
    Oprod2 <- myprod(delta, Omega2)
    Qprod1 <- myprod(delta, Q1)
    Qprod2 <- myprod(delta, Q2)

    if (het.reg) xi <- Qprod2/Qprod1
    phi1 <- sqrt(sigma1)
    phi2 <- sqrt(sigma2)
    p0 <- pargmaxV(0, phi1 = phi1, phi2 = phi2, xi = xi)
    if(is.nan(p0) || p0 < a2 || p0 > (1-a2)) {
      warning(sprintf("Confidence interval %d cannot be computed: P(argmax V <= 0) = %g", as.integer(i-1), round(p0, digits = 4)))
      upper[i-1] =  lower[i-1] = NA
    } else {
      
      ub <- lb <- 0
      while(pargmaxV(ub, phi1 = phi1, phi2 = phi2, xi = xi) < (1 - a2)) ub <- ub + 1000
      while(pargmaxV(lb, phi1 = phi1, phi2 = phi2, xi = xi) > a2) lb <- lb - 1000
      
      upper[i-1] <- uniroot(myfun, c(0, ub), level = (1-a2), xi = xi, phi1 = phi1, phi2 = phi2)$root
      lower[i-1] <- uniroot(myfun, c(lb, 0), level = a2, xi = xi, phi1 = phi1, phi2 = phi2)$root
      
      upper[i-1] <- upper[i-1] * phi1^2 / Qprod1
      lower[i-1] <- lower[i-1] * phi1^2 / Qprod1
    }#else
  }

  bp <- bp[-c(1, nbp + 2)]
  bp <- cbind(bp - ceiling(upper), bp, bp - floor(lower))
  a2 <- round(a2 * 100, digits = 1)
  colnames(bp) <- c(paste(a2, "%"), "breakpoints", paste(100 - a2, "%"))
  rownames(bp) <- 1:nbp

  return(bp)
}#computeConfInt


