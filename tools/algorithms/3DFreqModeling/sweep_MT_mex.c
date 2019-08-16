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
#define W prhs[4]
#define B prhs[3]
#define ROW_NORMS prhs[5]
#define OUTPUT plhs[0]
#define IDX1D(a,b,m,n) (a + b*m);
#define ATIMESB_RE(ar,ai,br,bi) ((ar)*(br) - (ai)*(bi)); 
#define ATIMESB_IM(ar,ai,br,bi) ((ar)*(bi) + (ai)*(br));

/*
 *
 * Matrix-vector product multi-threaded over right hand sides
 *
 * The N x N matrix A is given in diagonal storage format and the output is 
 *    y = A * x, where x is a N x P matrix.
 * 
 * 
 * Usage:
 *  y = sweep_MT_mex(coef,offset,x,b,w,row_norms_sq);
 *
 * Input:
 *   coef   - ndiags x N matrix of matrix coefficients of A
 *   offset - constant offset vector
 *   x      - N x P matrix
 *
 * Output:
 *   y      - multi-threaded carp sweep applied to x, with parameter w
 *  
 * Curt Da Silva, 2015, based on code by Mathias Louboutin and Art Petrenko
 *
 * To compile, run:
 *   mex -largeArrayDims sweep_MT_mex.c -DDEFINEUNIX -lmwblas CFLAGS="\$CFLAGS -fopenmp -O2" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 *
 *
 */


void init_zero(double *x, int n){
    int i;
    
    for(i=0;i<n;i++){
        x[i] = 0;
    }
}

void do_sweep(mwSize N, mwSize ndiags, mwSize P, double *coefr, double *coefi, double *offset, double *yr, double *yi, double* br, double* bi, double w, double * row_norms_sq ,mwSize nthreads) {      	 
    int m,i,j,k,kx,kout,xoffset,s;
    double atx_r,atx_i,c_r,c_i;
    
#pragma omp parallel for schedule(static,1) private(m,i,j,k,kx,kout,xoffset,atx_r,atx_i,c_r,c_i) num_threads(nthreads)
    for(s=0 ; s<P; s++)
	{       
	for(m=0; m<2*N ;m++ )
        {            
	    if (m >= N) { i = 2*N-m-1; } /* reverse sweep */
	    else { i = m; } /* forward sweep */

            kout = IDX1D(i,s,N,P);
	    atx_r = 0; atx_i = 0;
	   
	    /* atx = A(i,:) * y */
            for(j=0; j<ndiags; j++) {
                xoffset = offset[j];                
                if(i+xoffset >= 0 && i+xoffset < N){
                    k = IDX1D(j,i,ndiags,N);                        
		    kx = IDX1D(i+xoffset,s,N,P); 
                    atx_r += ATIMESB_RE(coefr[k],coefi[k],yr[kx],yi[kx]);
                    atx_i += ATIMESB_IM(coefr[k],coefi[k],yr[kx],yi[kx]);
                }                                
            }	    
	    
	    /* c = w*(b(i) - A(i,:)* y)/|A(i,:)|^2 */
	    c_r = w*(br[kout] - atx_r)/row_norms_sq[i];
	    c_i = w*(bi[kout] - atx_i)/row_norms_sq[i];
	    
	    /* y <- y + c * conj(A(i,:)) */
	    for(j=0; j<ndiags; j++) {
		xoffset = offset[j];                
                if(i+xoffset >= 0 && i+xoffset < N){
		    k = IDX1D(j,i,ndiags,N);                        
		    kx = IDX1D(i+xoffset,s,N,P);
		    yr[kx] += ATIMESB_RE(coefr[k],-coefi[k],c_r,c_i);
		    yi[kx] += ATIMESB_IM(coefr[k],-coefi[k],c_r,c_i);
		}
	    }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    mwSize N, P,ndiags;        
    mwSize n_threads = 1;
    
    
    double *coefr, *coefi, *offset, *yr, *yi, *xr, *xi,*br,*bi,*row_norms_sq;     
    double w = 1.5;
    int Si_alloc=0, xi_alloc=0,bi_alloc=0;        
    
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
   
    n_threads = P;
    xr  = mxGetPr(X);
    if(mxIsComplex(X)){
        xi  = mxGetPi(X);
    }
    else{
        xi = mxCalloc(N*P,sizeof(double));
        init_zero(xi,N*P);
        xi_alloc = 1;
    }

    br = mxGetPr(B);
    if(mxIsComplex(B)){
	bi = mxGetPi(B);
    }
    else{
	bi = mxCalloc(N*P,sizeof(double));
	init_zero(bi,N*P);
	bi_alloc = 1;
    }
    row_norms_sq = mxGetPr(ROW_NORMS);
    if(mxGetM(ROW_NORMS) != N)
    {
	mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:HMvp_MT_mex:NumElements",
			  "The size of the row norms vector and the rows size of the matrix A must match");
    }
    
    w = *(mxGetPr(W));
    
    /* define output vector y and initialize with input vector x.*/
    OUTPUT   = mxCreateDoubleMatrix(N, P, mxCOMPLEX);    
    yr       = mxGetPr(OUTPUT);
    yi       = mxGetPi(OUTPUT);

    memcpy((void*)yr,(void*)xr,sizeof(double)*N*P);
    memcpy((void*)yi,(void*)xi,sizeof(double)*N*P);
    
    
    /*multi-threaded sweep*/
    do_sweep(N, ndiags,  P, coefr,coefi,offset,yr,yi,br,bi,w,row_norms_sq,n_threads);    
    
    if (Si_alloc){
        mxFree(coefi);
    }
    if (xi_alloc){
        mxFree(xi);
    }
    if (bi_alloc){
	mxFree(bi);
    }
    
    return;
}
