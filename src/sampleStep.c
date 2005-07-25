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
  int i,j,imax;

  for(i=0; i<n; i++) 
    LOGICAL(ans)[i] = TRUE;

  for(i=1; i<n; i++) {
    if(x[i-1]>x[i])
      error("Elements of x must be in ascending order.");
  }

  imax = n-1;
  for(i=0; i<imax; ) {
    for(j=i+1; (j<n) && (x[j]-x[i]<step); j++)
      LOGICAL(ans)[j] = FALSE;
    i=j;  
  }
}
/*-----------------------------------------------------------------
sampleSteps
------------------------------------------------------------------*/
SEXP sampleStep(SEXP ax, SEXP astep) 
{
  SEXP ans;   /* return value    */
  double *x;
  double step;
  int n;      /* length of x */

  /* check input arguments */
  if((!isReal(ax)))
    error("Invalid argument 'ax', must be real."); 
  x    = REAL(ax);
  n    = length(ax);

  if(!isReal(astep) || length(astep)!=1)
      error("Invalid argument 'astep', must be real of length 1.");
  step = REAL(astep)[0];

  /* return value: a logical vector of length n */
  PROTECT(ans = allocVector(LGLSXP, n));

  sampleStep_c(x, n, step, ans);
  
  UNPROTECT(1); /* done with ans */
  return ans;
}

