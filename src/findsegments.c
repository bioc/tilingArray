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
#define MAT_ELT(x, i, j, nrow) x[i+(j)*(nrow)]

/*  Global variables */
double *G;  /* cost matrix                                        */
int maxk;   /* number of rows of G:    maximum length of segments */
int n;      /* number of columns of G: number of data points      */

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
  accessor function for matrix G
  G[k, i] = G[k+i*maxk] is the cost of segment from i to i+k,
    including these endpoints
  fG(i, j) calculates the cost for segment from i to j-1, excluding
    j. (j=i+k+1).
-----------------------------------------------------------------*/
static R_INLINE double fG(int i, int j) {
    int k;
    /* Rprintf("fG: %4d %4d\n", i, j); */
#ifdef DEBUG
    if((i<0) || (i>=n) || (j<0) || (j>=n))
	error("Illegal value for i or j.");
#endif
    k = j-i-1;
    return ((k>=0) && (k<maxk)) ? MAT_ELT(G, k, i, maxk) : R_PosInf;
}


/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.  This is the workhorse routine with C interface.
   Note that all array indices here start at 0 and run to size(array)-1.
   At the end we add 1 to the result indices in matrix 'th'
-----------------------------------------------------------------*/
void findsegments_dp(double* J, int* th, int maxcp) {
    int i, imin, cp, j;
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

    /* initialize for cp=0: mI[0, i] is simply fG(0, i) */
    for(i=0; i<n; i++)
	MAT_ELT(mI, 0, i, maxcp) = fG(0, i);

    for (cp=1; cp<maxcp; cp++) {	
      /*  Best segmentation with cp change points from 0 to j 
      is found from considering best segmentations from 0 to i (<j)
      with cp-1 segments, plus cost of segment from i to j */
      for (j=0; j<n; j++) {   
	  zmin = R_PosInf;
	  imin = j;
          /* find the best change point between 0 and j-1 */ 
	  for (i=1; i<j; i++) { 
	      z = MAT_ELT(mI, cp-1, i, maxcp) + fG(i, j);
              /* Rprintf("%2d %2d %2d %6g %6g %6g\n", 
		 cp, j, i, mI[cp-1 + i*maxcp], fG(i, j), z); */
	      if(z<zmin) {
		  zmin = z;
		  imin = i;
              }
	  } /* for i */	  
	  MAT_ELT(mI, cp,   j, maxcp  ) = zmin;
	  MAT_ELT(mt, cp-1, j, maxcp-1) = imin;
      } /* for j */
    } /* for cp */

    /* print_matrix_double(mI, maxcp, n, "mI");
       print_matrix_int(mt, maxcp-1, n, "mt"); */
   
    /* th: elements 0...cp-1 of the cp-th row of matrix th contain
       the cp change points; element cp has value n, which corresponds
       to a changepoint at the rightmost point */
    for(cp=0;  cp<maxcp; cp++) {
	/* Calculate J, the log-likelihood. TO DO: constant factors, sqrt(2*pi) */
	J[cp] = -log(MAT_ELT(mI, cp, n-1, maxcp) / n);

	for(j=cp+1; j<maxcp; j++)
	    MAT_ELT(th, cp, j, maxcp) = -1;

	/* Backtrack to get th */
        /* In the following loop, i is always the changepoint to the right */
        MAT_ELT(th, cp, cp, maxcp) = i = n;  /* note the chained assignment */
	for(j=cp-1; j>=0; j--) {
	    Rprintf("cp=%4d j=%4d i=%4d\n", cp, j, i); 
#ifdef DEBUG
	    if((i<1)||(i>n)) {
	       error("Illegal value for i.");
	    }
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
SEXP findsegments(SEXP _G, SEXP _maxcp) 
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
  
  /* J */
  PROTECT(J   = allocVector(REALSXP, maxcp));

  /* th */
  PROTECT(th    = allocVector(INTSXP, maxcp*maxcp));
  PROTECT(dimth = allocVector(INTSXP, 2));
  INTEGER(dimth)[0] = INTEGER(dimth)[1] = maxcp;
  setAttrib(th, R_DimSymbol, dimth);

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
    { "findsegments", ( DL_FUNC ) &findsegments, 2,
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

