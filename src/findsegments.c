/*
 * Copyright W. Huber 2005, all rights reserved
 */
 
/* 
 * Most of this code was written on Detroit airport waiting for 
 * a delayed Northwest flight connection to Columbus.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h>

#define DEBUG
#undef VERBOSE

#define MAT_ELT(x, i, j, nrow) x[i+(j)*(nrow)]

/*  Global variables */
double *G;  /* cost matrix                                        */
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

/*-----------------------------------------------------------------
  Accessor function for matrix G
  G[k, i] is the cost of segment from i to i+k, including these endpoints
  fG(i, j) calculates the cost for segment from i to j-1, excluding
    j. (j=i+k+1).
-----------------------------------------------------------------*/
static R_INLINE double fG(int i, int j) {
    int k;
    if(i==j)
      return((double)0);
    k = j-i-1;
#ifdef DEBUG
    if((i<0) || (i>=n) || (k<0) || (k>=maxk)) {
      Rprintf("i=%d j=%d k=%d\n", i, j, k);
      error("Illegal value.");
    }
#endif
    return(MAT_ELT(G, k, i, maxk));
}


/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.  This is the workhorse routine with C interface.
   Note that all array indices here start at 0 and run to size(array)-1.
   At the end we add 1 to the result indices in matrix 'th'
-----------------------------------------------------------------*/
void findsegments_dp(double* J, int* th, int maxcp) {
    int i, imin, i0, cp, j;
    double z, zmin;
    double *mI;
    int * mt;
    
    /* mI[cp, i] is the optimal cost of segmentation from 0 to i 
       with cp change points */
    /* mt[cp-1, i] is the index (0...n-1) of the rightmost changepoint 
       in the optimal segmentation from 0 to i with cp change points;
       the whole segmentation can then be reconstructed from recursing 
       through this matrix */
    mI = (double*) R_alloc(  maxcp  *n, sizeof(double));
    mt = (int*)    R_alloc((maxcp-1)*n, sizeof(int));

    /* initialize for cp=0: mI[0, j] is simply fG(0, j) */
    for(j=0; j<maxk; j++)
	MAT_ELT(mI, 0, j, maxcp) = fG(0, j);
    for(j=maxk; j<n; j++)
	MAT_ELT(mI, 0, j, maxcp) = R_PosInf;

    for (cp=1; cp<maxcp; cp++) {	
      if(verbose>=2)
	Rprintf("findsegments_dp: cp=%d.\n", cp);
      /*  Best segmentation with cp change points from 0 to j 
      is found from considering best segmentations from 0 to i (<j)
      with cp-1 segments, plus cost of segment from i to j. 
      And j-maxk <= i <= j */
      for (j=0; j<n; j++) {   
	  zmin = R_PosInf;
	  imin = j;
          /* find the best change point between 0 and j-1 */ 
          i0 = j-maxk;
          if(i0<0) i0=0;
	  for (i=i0; i<=j; i++) { 
#ifdef VERBOSE              
	      Rprintf("%2d %2d %2d %6g %6g", 
		      cp, j, i, MAT_ELT(mI, cp-1, i, maxcp), fG(i, j));
#endif
	      z = MAT_ELT(mI, cp-1, i, maxcp);
              if (finite(z))
		  z += fG(i, j);
#ifdef VERBOSE              
	      Rprintf(" %6g\n", z);
#endif
	     if(z<zmin) {
		  zmin = z;
		  imin = i;
	     } /* if z */
	  } /* for i */	  
	  MAT_ELT(mI, cp,   j, maxcp  ) = zmin;
	  MAT_ELT(mt, cp-1, j, maxcp-1) = imin;
      } /* for j */
    } /* for cp */

#ifdef VERBOSE              
       print_matrix_double(mI, maxcp, n, "mI");
       print_matrix_int(mt, maxcp-1, n, "mt"); 
#endif
   
    /* th: elements 0...cp-1 of the cp-th row of matrix th contain
       the cp change points; element cp has value n, which corresponds
       to a changepoint at the rightmost point */
    for(cp=0;  cp<maxcp; cp++) {
        /* Calculate J, the log-likelihood. TO DO: constant factors, sqrt(2*pi) */
        z = MAT_ELT(mI, cp, n-1, maxcp);
	J[cp] = finite(z) ? -log(z/n) : R_NegInf;

	for(j=cp+1; j<maxcp; j++)
	    MAT_ELT(th, cp, j, maxcp) = -1;

	/* Backtrack to get th */
        /* In the following loop, i is always the changepoint to the right */
        MAT_ELT(th, cp, cp, maxcp) = i = n;  /* note the chained assignment */
	for(j=cp-1; j>=0; j--) {
#ifdef VERBOSE              
	    Rprintf("cp=%4d j=%4d i=%4d\n", cp, j, i); 
#endif
#ifdef DEBUG
	    if((i<1)||(i>n))
	       error("Illegal value for i.");
#endif
            /* note the chained assignment */
	    MAT_ELT(th, cp, j, maxcp) = i = MAT_ELT(mt, j, i-1, maxcp-1);
	}
    }

    /* add 1 to all elements of th since in R array indices start at 1,
       while here they were from 0 */
    for(i=0; i<cp*cp; i++) 
       th[i] += 1; 

    return;
} 

/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.
   This is the R interface:
   G : numeric matrix (cost matrix)
   maxcp: integer scalar (maximum number of segments
------------------------------------------------------------------*/
SEXP findsegments(SEXP _G, SEXP _maxcp, SEXP _verbose) 
{
  SEXP dimG;  /* dimensions of G */
  SEXP res;   /* return value    */
  SEXP J, th, dimth, namesres;  /* for the return value */
  int maxcp;

  /* check input arguments */
  PROTECT(dimG = getAttrib(_G, R_DimSymbol));
 
  if((!isReal(_G)) | isNull(dimG) | (LENGTH(dimG)!=2))
    error("Invalid argument '_G', must be a real matrix."); 
  G    = REAL(_G);
  maxk = INTEGER(dimG)[0];
  n    = INTEGER(dimG)[1];
  UNPROTECT(1);  /* done with dimG */

  if(!isInteger(_maxcp) | length(_maxcp)!=1)
      error("'_maxcp' must be integer of length 1.");
  maxcp = INTEGER(_maxcp)[0];
  
  if(!isInteger(_verbose) | length(_verbose)!=1)
      error("'_verbose' must be integer of length 1.");
  verbose = INTEGER(_verbose)[0];

  /* J */
  PROTECT(J   = allocVector(REALSXP, maxcp));

  /* th */
  PROTECT(th    = allocVector(INTSXP, maxcp*maxcp));
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
  return(res);
}


extern void R_init_tilingArray( DllInfo *info );
extern void R_unload_tilingArray( DllInfo *info );

/* Registration information for DLL */
static R_CallMethodDef callMethods[] = {
    { "findsegments", ( DL_FUNC ) &findsegments, 3,
        /* { REALSXP, INTSXP, INTSXP } */ },
    { NULL, NULL, 0 }
};

void R_init_tilingArray( DllInfo *info ) {
    /* R_useDynamicSymbols( dll, FALSE );  don't know what this command does */
    R_registerRoutines( info, NULL, callMethods, NULL, NULL );
}

void R_unload_tilingArray( DllInfo *info ) {
  /* Release resources. */
}

