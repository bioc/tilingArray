/*
 * Copyright W. Huber 2005-2009, all rights reserved
 */
 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>

#define DEBUG

#define MAT_ELT(x, i, j, nrow) x[(long)(i)+(long)(j)*(long)(nrow)]

/* from sampleSteps.c */
SEXP sampleStep(SEXP ax, SEXP astep);

/*  Global variables */
double *G;  /* cost matrix */
int maxk;   /* number of rows of G:    maximum length of segments */
int n;      /* number of columns of G: number of data points      */
int verbose;

/*--------------------------------------------------
  For debugging
---------------------------------------------------*/
void print_matrix_double(double* x, int nrow, int ncol, char *s) {
    int i, j;
    Rprintf("%s:\n", s);
    for(i=0; i<nrow; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<ncol; j++) 
	    Rprintf("%7g ", MAT_ELT(x, i, j, nrow));
        Rprintf("\n");
    }
}

void print_matrix_int(int* x, int nrow, int ncol, char *s) {
    int i, j;
    Rprintf("%s:\n", s);
    for(i=0; i<nrow; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<ncol; j++) 
	    Rprintf("%3d ", MAT_ELT(x, i, j, nrow));
        Rprintf("\n");
    }
}

R_len_t matrix_length(int i, int j){
  double x;
  if ((i<0) || (j<0))
    error("Negative values not allowed in 'matrix_length'.");
  x = (double) i * (double) j;
  if(x >  R_LEN_T_MAX)
    error("Please do not try to allocate a vector whose length is greater than R_LEN_T_MAX.");
  return((R_len_t) x);
}

/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.  This is the workhorse routine with C interface.
   Note that all array indices here start at 0 and run to size(array)-1.
   At the end we add 1 to the result indices in matrix 'th'
-----------------------------------------------------------------*/
void findsegments_dp(double* J, int* th, int maxcp) {
    int i, imin, cp, j, k, k0;
    double z, zmin;
    double *mI;
    int * mt;
    SEXP v1, v2;    

    if(verbose>=2)
	Rprintf("In findsegments_dp: cp=      ");

    /* G[k, i] is the cost of segment from i to i+k, including these 
       endpoints */
    /* mI[i, cp] is the optimal cost of segmentation from 0 to i 
       with cp change points */
    /* mt[i, cp-1] is the index (0...n-1) of the rightmost changepoint 
       in the optimal segmentation from 0 to i with cp change points;
       the whole segmentation can then be reconstructed from recursing 
       through this matrix */

    PROTECT(v1 = allocVector(REALSXP, matrix_length(maxcp, n))); 
    mI = REAL(v1);

    PROTECT(v2 = allocVector(INTSXP, matrix_length(maxcp-1, n))); 
    mt = INTEGER(v2);

    /* initialize for cp=0: mI[k, 0] is simply G[k, 0] */
    for(k=0; k<maxk; k++)
	mI[k] = G[k];
    for(k=maxk; k<n; k++)
	mI[k] = R_PosInf;

    for (cp=1; cp<maxcp; cp++) {	
      if(verbose>=2)
	Rprintf("\b\b\b\b\b\b%6d", cp);
      /*  Best segmentation with cp change points from 0 to j 
      is found from considering best segmentations from 0 to j-k-1
      with cp-1 segments, plus cost of segment from j-k to j. */
      for (j=0; j<n; j++) {   
	  zmin = R_PosInf;
	  imin = j;
          /* find the best change point between 0 and j-1 */ 
          k0 = (j<maxk) ? j : maxk;
	  for (k=0; k<k0; k++) { 
              /* Best segmentation from 0 to j-k-1 */
	      z = MAT_ELT(mI, j-k-1, cp-1, n);
              if (R_FINITE(z))
		  z += MAT_ELT(G, k, j-k, maxk);
                  /* Cost of segment from j-k to j */
	      if(z<zmin) {
		  zmin = z;
		  imin = j-k;
	     } /* if z */
	  } /* for k */	  
	  MAT_ELT(mI, j,   cp, n) = zmin;
	  MAT_ELT(mt, j, cp-1, n) = imin;
      } /* for j */
      /* to allow user interrupts */
      R_CheckUserInterrupt();
    } /* for cp */

    if(verbose>=2)
	Rprintf("\n");

    /* th: elements 0...cp-1 of the cp-th row of matrix th contain
       the cp change points; element cp has value n, which corresponds
       to a changepoint at the rightmost point */
    for(cp=0;  cp<maxcp; cp++) {
        /* Calculate J, the log-likelihood. */
        z = MAT_ELT(mI, n-1, cp, n);
	J[cp] = R_FINITE(z) ? -(double)n/2.0*(1+log(2*M_PI)+log(z/n)) : R_NegInf;

	for(j=cp+1; j<maxcp; j++)
	    MAT_ELT(th, cp, j, maxcp) = -1;

	/* Backtrack to get th */
        /* In the following loop, i is always the changepoint to the right */
        MAT_ELT(th, cp, cp, maxcp) = i = n;  /* note the chained assignment */
	for(j=cp-1; j>=0; j--) {
#ifdef DEBUG
	    if((i<1)||(i>n))
	       error("Illegal value for i.");
#endif
            /* note the chained assignment */
	    MAT_ELT(th, cp, j, maxcp) = i = MAT_ELT(mt, i-1, j, n);
	}
    }

    UNPROTECT(2);

    /* add 1 to all elements of th since in R array indices start at 1,
       while here they were from 0 */
    for(cp=0; cp<maxcp*maxcp; cp++) 
       th[cp] += 1; 

} 

/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.
   This is the R interface:
   G : numeric matrix (cost matrix)
   maxcp: integer scalar (maximum number of segments
------------------------------------------------------------------*/
SEXP findsegments(SEXP aG, SEXP amaxcp, SEXP averbose) 
{
  SEXP dimG;  /* dimensions of G */
  SEXP res;   /* return value    */
  SEXP J, th, dimth, namesres;  /* for the return value */
  int maxcp;

  /* check input arguments */
  PROTECT(dimG = getAttrib(aG, R_DimSymbol));
 
  if((!isReal(aG)) || isNull(dimG) || (LENGTH(dimG)!=2))
    error("Invalid argument 'aG', must be a real matrix."); 
  G    = REAL(aG);
  maxk = INTEGER(dimG)[0];
  n    = INTEGER(dimG)[1];
  UNPROTECT(1);  /* done with dimG */

  if(!isInteger(amaxcp) || length(amaxcp)!=1)
      error("'amaxcp' must be integer of length 1.");
  maxcp = INTEGER(amaxcp)[0];
  
  if(!isInteger(averbose) || length(averbose)!=1)
      error("'averbose' must be integer of length 1.");
  verbose = INTEGER(averbose)[0];

  /* J */
  PROTECT(J   = allocVector(REALSXP, maxcp));

  /* th */
  PROTECT(th    = allocVector(INTSXP, matrix_length(maxcp, maxcp)));
  PROTECT(dimth = allocVector(INTSXP, 2));
  INTEGER(dimth)[0] = INTEGER(dimth)[1] = maxcp;
  setAttrib(th, R_DimSymbol, dimth);

  if(verbose>=2)
    Rprintf("In C code now, maxk=%d, n=%d, maxcp=%d\n", maxk, n, maxcp);

  findsegments_dp(REAL(J), INTEGER(th), maxcp);

  /* return value: a list with two elements, th and J */
  PROTECT(res = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, J);
  SET_VECTOR_ELT(res, 1, th);

  PROTECT(namesres = allocVector(STRSXP, 2));
  SET_STRING_ELT(namesres, 0, mkChar("J"));
  SET_STRING_ELT(namesres, 1, mkChar("th"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(5); /* done with res, namesres, J, th, dimth */
  return res;
}


extern void R_init_tilingArray( DllInfo *info );
extern void R_unload_tilingArray( DllInfo *info );

/* Registration information for DLL */
static R_CallMethodDef callMethods[] = {
    { "findsegments", ( DL_FUNC ) &findsegments, 3,
        /* { REALSXP, INTSXP, INTSXP } */ },
    { NULL, NULL, 0 },
    { "sampleStep", ( DL_FUNC ) &sampleStep, 2,
        /* { REALSXP, REALSXP } */ },
    { NULL, NULL, 0 }
};

void R_init_tilingArray( DllInfo *info ) {
    /* R_useDynamicSymbols( dll, FALSE );  don't know what this command does */
    R_registerRoutines( info, NULL, callMethods, NULL, NULL );
}

void R_unload_tilingArray( DllInfo *info ) {
  /* Release resources. */
}

