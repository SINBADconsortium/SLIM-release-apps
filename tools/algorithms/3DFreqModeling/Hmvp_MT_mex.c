#include <stdlib.h> /* for getenv */
#include <stddef.h> /* for size_t type */
#include <string.h> /* for memcpy */
#include <pthread.h> /* for threading */
#include <omp.h>

#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <assert.h>

#define COEF prhs[0]
#define OFFSET prhs[1]
#define X prhs[2]
#define NTHREADS prhs[3]
#define OUTPUT plhs[0]
#define IDX1D(a,b,m,n) (a + b*m);
#define ATB_RE(ar,ai,br,bi) (ar*br - ai*bi); 
#define ATB_IM(ar,ai,br,bi) (ar*bi + ai*br);

/*
 *
 * Matrix-vector product multi-threaded over right hand sides
 *
 * The N x N matrix A is given in diagonal storage format and the output is 
 *    y = A * x, where x is a N x P matrix. Uses P threads.
 * 
 * 
 * Usage:
 *  y = Hmvp_MT_mex(coef,offset,x);
 *
 * Input:
 *   coef   - ndiags x N matrix of matrix coefficients of A
 *   offset - constant offset vector
 *   x      - N x P matrix
 *
 * Output:
 *   y      - A * x
 *  
 * Curt Da Silva, 2015
 *
 * To compile, run:
 *   mex -O -largeArrayDims Hmvp_MT_mex.c -DDEFINEUNIX -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 *
 */


void init_zero(double *x, int n){
    int i;
    
    for(i=0;i<n;i++){
        x[i] = 0;
    }
}


void do_Hmvp(mwSize N, mwSize ndiags, mwSize P, double *coefr, double *coefi, double *offset, double *yr, double *yi, double *xr, double *xi,mwSize nthreads) {      	 
    int i,j,k,kx,kout,xoffset,s,t;
    
#pragma omp parallel for schedule(static) private(s,i,j,k,kx,kout,xoffset) num_threads(nthreads)
    for(t=0 ; t<P*N; t++)
    {       
	s = t/N; i = t%N;
	
	kout = IDX1D(i,s,N,P);
	for(j=0; j<ndiags; j++) {
	    xoffset = i+offset[j];                
	    if(xoffset >= 0 && xoffset < N){
		k = IDX1D(j,i,ndiags,N);
		kx = IDX1D(xoffset,s,N,P);                     
		yr[kout] += ATB_RE(coefr[k],coefi[k],xr[kx],xi[kx]);
		yi[kout] += ATB_IM(coefr[k],coefi[k],xr[kx],xi[kx]);
	    }                                
	}
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    mwSize N, P,ndiags;        
    mwSize n_threads;
    
    
    double *coefr, *coefi, *offset, *yr, *yi, *xr, *xi;     
    int Si_alloc=0, xi_alloc=0;
        
    char *n_threads_str = NULL;
    
    /* read input, initialize complex part to zero if input is real.*/
    N        = mxGetN(COEF);
    ndiags   = mxGetM(COEF);
    coefr    = mxGetPr(COEF);
    if(mxIsComplex(COEF)){
        coefi  = mxGetPi(COEF);
    }
    else{
        coefi = mxCalloc(N*ndiags,sizeof(double));
        init_zero(coefi,N*ndiags);
        Si_alloc = 1;
    }
    offset      = mxGetPr(OFFSET);    
    P           = mxGetN(X);
    
    if(mxGetM(X) != N)
    {
        mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:HMvp_MT_mex:NumElements",
                "The row size of the vector x and the column size of the matrix A must match");
    }
                   
    n_threads = mxGetScalar(NTHREADS);
    xr  = mxGetPr(X);
    if(mxIsComplex(X)){
        xi  = mxGetPi(X);
    }
    else{
        xi = mxCalloc(N*P,sizeof(double));
        init_zero(xi,N*P);
        xi_alloc = 1;
    }
    
    /* define output vector y and initialize with input vector x.*/
    OUTPUT   = mxCreateDoubleMatrix(N, P, mxCOMPLEX);    
    yr       = mxGetPr(OUTPUT);
    yi       = mxGetPi(OUTPUT);
    init_zero(yr,N*P);
    init_zero(yi,N*P);
    
    /*Hmvp*/
    do_Hmvp(N, ndiags,  P, coefr,coefi,offset,yr,yi,xr,xi,n_threads);    
    
    if (Si_alloc){
        mxFree(coefi);
    }
    if (xi_alloc){
        mxFree(xi);
    }
    
    return;
}
