/*
 * Copyright W. Huber 2005, all rights reserved
 */
 
/* 
 * Most of this code was written on Detroit airport waiting for 
 * a delayed Northwest flight connection to Columbus.
 */

#include <R.h>
#include <Rinternals.h>
/* #include "R_ext/Applic.h"  */
#include <R_ext/Rdynload.h>

/* #include <string.h> */
#include <stdlib.h>

#define DEBUG

static R_INLINE double sqr(double x) { 
  return(x*x); 
}

/*  Global variables */
int maxk; /* number of rows of G :    maximum length of segments */
int n;    /* number of columns of G : number of data points      */
double *G;   /* cost matrix */

/*--------------------------------------------------
  For debugging
---------------------------------------------------*/
void print_matrix_double(double* x, int n1, int n2, char *s) {
    int i, j;
    Rprintf("%s:\n", s);
    for(i=0; i<n1; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<n2; j++) 
	    Rprintf("%7g ", x[i+j*n1]);
        Rprintf("\n");
    }
}
void print_matrix_int(int* x, int n1, int n2, char *s) {
    int i, j;
    Rprintf("%s:\n", s);
    for(i=0; i<n1; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<n2; j++) 
	    Rprintf("%3d ", x[i+j*n1]);
        Rprintf("\n");
    }
}

/*-----------------------------------------------------------------
  accessor function for matrix G
  G[k, i] = G[k+i*maxk] is the cost of segment from i to i+k,
    including these endpoints
  fG(i, j) calculates the cost for segment from i to j-1, including
    these endpoints (j=i+k+1).
-----------------------------------------------------------------*/
static R_INLINE double fG(int i, int j) {
    int k;
    /* Rprintf("fG: %4d %4d\n", i, j); */
#ifdef DEBUG
    if((i<0) || (i>=n) || (j<0) || (j>=n))
	error("Illegal value for i or j.");
#endif
    k = j-i-1;
    return ((k>=0) && (k<maxk)) ? G[k+i*maxk] : R_PosInf;
}


/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.  This is the workhorse routine with C interface
-----------------------------------------------------------------*/
void findsegments_dp(double* J, int* th, int Km) {
    int i, imin, k, j;
    double z, zmin;
    double *mI;
    int * mt;
    
    /* mI[k, i] is the optimal cost of segmentation from 0 to i 
       with k segments */
    mI = (double*) R_alloc(  Km  *n, sizeof(double));
    mt = (int*)    R_alloc((Km-1)*n, sizeof(int));

    /* initialize for k=0: mI[0, i] is simply fG(0, i) */
    for(i=0; i<n; i++)
	mI[i*Km] = fG(0, i);

    /* k: the number of change points */
    for (k=1; k<Km; k++) {	
      /*  Best segmentation with k segments from 0 to j 
      is found from considering best segmentations from 0 to i (<j)
      with k-1 segments, plus cost of segment from i to j */
      for (j=0; j<n; j++) {   
	  zmin = R_PosInf;
	  imin = -1;
          /* find the best jump point between 0 and j-1 */ 
	  for (i=0; i<j; i++) { 
	      z = mI[k-1 + i*Km] + fG(i, j);
              /* Rprintf("%2d %2d %2d %6g %6g %6g\n", 
		 k, j, i, mI[k-1 + i*Km], fG(i, j), z); */
	      if(z<zmin) {
		  zmin = z;
		  imin = i;
              }
	  } /* for i */	  
	  mI[k   + j*Km    ] = zmin;
	  mt[k-1 + j*(Km-1)] = imin;
      } /* for j */
    } /* for k */

    print_matrix_double(mI, Km, n, "mI");
    print_matrix_int(mt, Km-1, n, "mt");
   
    /* elements 0...k-1 of the k-th row of matrix th contain
       the k change points; element k has value n, which corresponds
       to a changepoint at the rightmost point */
    for(k=0;  k<Km; k++) {
	/* Calculate J */
	J[k] = log(mI[k + (n-1)*Km]/n);
	/* Backtrack to get th */
        th[k + k*Km] = n;
	for(j=k+1; j<Km; j++)
	    th[k + j*Km] = -1;
	for(j=k-1; j>=0; j--) {
	    /* i  is the changepoint to the right */
            i = th[k+(j+1)*Km]-1;
            Rprintf("%2d %2d %2d %3d\n", k, j, i, j + i*(Km-1));
#ifdef DEBUG
	    if((i<0)||(i>=n))
		error("Illegal value for i.");
#endif
	    th[k + j*Km] = mt[j + i*(Km-1)];
	}
    }

    return;
} 

/*-----------------------------------------------------------------
   Find segments using the dynamic programming algorithm of Picard
   et al.
   This is the R interface:
   G : numeric matrix (cost matrix)
   Km: integer scalar (maximum number of segments
------------------------------------------------------------------*/
SEXP findsegments(SEXP _G, SEXP _Km) 
{
  SEXP dimG;  /* dimensions of G */
  SEXP res;   /* return value    */
  SEXP J, th, dimth, namesres;  /* for the return value */
  int Km;

  /* check input arguments */
  PROTECT(dimG = getAttrib(_G, R_DimSymbol));
 
  if((!isReal(_G)) | isNull(dimG) | (LENGTH(dimG)!=2))
    error("Invalid argument '_G', must be a real matrix."); 
  G    = REAL(_G);
  maxk = INTEGER(dimG)[0];
  n    = INTEGER(dimG)[1];
  UNPROTECT(1);  /* done with dimG */

  if(!isInteger(_Km) | length(_Km)!=1)
      error("'_Km' must be integer of length 1.");
  Km = INTEGER(_Km)[0];
  
  /* J */
  PROTECT(J   = allocVector(REALSXP, Km));

  /* th */
  PROTECT(th    = allocVector(INTSXP, Km*Km));
  PROTECT(dimth = allocVector(INTSXP, 2));
  INTEGER(dimth)[0] = INTEGER(dimth)[1] = Km;
  setAttrib(th, R_DimSymbol, dimth);

  findsegments_dp(REAL(J), INTEGER(th), Km);

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

