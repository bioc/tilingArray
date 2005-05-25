/*
 * Copyright W. Huber 2005, all rights reserved
 */
 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>

void sampleStep_c(double* x, int n, double step, SEXP ans)
{
  int i,j,imax, jmax;

  for(i=0; i<n; i++) 
    LOGICAL(ans)[i] = TRUE;

  for(i=1; i<n; i++) {
    if(x[i-1]>x[i])
      error("Elements of x must be in ascending order.");
  }

  imax = n-2;
  jmax = n-1;
  for(i=0; i<imax; ) {
    for(j=i+1; (j<jmax) && (x[j+1]-x[i]<step); j++)
      LOGICAL(ans)[j] = FALSE;
    i=j;  
  }
}
/*-----------------------------------------------------------------
sampleSteps
------------------------------------------------------------------*/
SEXP sampleStep(SEXP _x, SEXP _step) 
{
  SEXP ans;   /* return value    */
  double *x;
  double step;
  int n;      /* length of x */

  /* check input arguments */
  if((!isReal(_x)))
    error("Invalid argument '_x', must be real."); 
  x    = REAL(_x);
  n    = length(_x);

  if(!isReal(_step) || length(_step)!=1)
      error("Invalid argument '_step', must be real of length 1.");
  step = REAL(_step)[0];

  /* return value: a logical vector of length n */
  PROTECT(ans = allocVector(LGLSXP, n));

  sampleStep_c(x, n, step, ans);
  
  UNPROTECT(1); /* done with ans */
  return ans;
}

