#include <stdlib.h> /* for getenv */
#include <stddef.h> /* for size_t type */
#include <string.h> /* for memcpy */
#include <omp.h>

#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "Helm3d_27pt.h"

/*
 *
 * 27pt stencil Helmholtz matrix-vector product.  
 * 
 * 
 * Usage:
 *  y = Helm3dmvp_mex(wn,h,n,nmpl,x,nthreads);
 *
 * Input:
 *   wn       - wavenumber evaluated on pml grid
 *   h        - [hx,hy,hz] spacings in each direction
 *   n        - [nx,ny,nz] points in each dimension, including pml
 *   npml     - 2 x 3 matrix denoting the number of pml points in each dimension                    
 *   x        - nx*ny*nz x 1 vector
 *   nthreads - number of threads
 *
 * Output:
 *   y      - H * x
 *  
 * Curt Da Silva, 2015
 *
 * To compile, run:
 *   mex -O -largeArrayDims Helm3dmvp_mex.c -DDEFINEUNIX -lmwblas CFLAGS="\$CFLAGS -fopenmp -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 * Compiler flags:
 *   -DADJ   - adjoint matrix-vector product
 *   -DDERIV - derivative matrix mvp
 */


#define WN prhs[0]
#define H prhs[1]
#define N prhs[2]
#define NPML prhs[3]
#define X prhs[4]
#define NTHREADS prhs[5]
#define OUTPUT plhs[0]

void init_zero(double *x, int n){
    int i;
    
    for(i=0;i<n;i++){
        x[i] = 0;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    mwSize n_threads = 1;
        
    double * wnr, *wni, *h, *n, *npml, *yr, *yi, *xr, *xi;     
    mwSize numel;
    char *n_threads_str = NULL;       
    int xi_alloc = 0;
    
    wnr   = mxGetPr(WN);
    wni   = mxGetPi(WN);
    h     = mxGetPr(H);
    n     = mxGetPr(N);
    npml  = mxGetPr(NPML);                       
    
    n_threads = mxGetScalar(NTHREADS);

    if(mxGetM(N) != 3 && mxGetN(N) != 3)
    {
	mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:Helm3dmvp_mex:Nsize",
			  "n must be a 3-length vector");
    }

    numel = n[0]*n[1]*n[2];
    if(mxGetM(X) != numel) {
        mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:Helm3dkaczswp_mex:Xsize",
			  "x must have nx*ny*nz elements");
    }
    xr  = mxGetPr(X);
    if(mxIsComplex(X)){
        xi  = mxGetPi(X);
    }
    else{
        xi = mxCalloc(numel,sizeof(double));
	xi_alloc = 1;
    }
       
     /* define output vector y and initialize with input vector x.*/
    OUTPUT   = mxCreateDoubleMatrix(numel, 1, mxCOMPLEX);    
    yr       = mxGetPr(OUTPUT);
    yi       = mxGetPi(OUTPUT);
    init_zero(yr,numel);
    init_zero(yi,numel);

    do_Hmvp( wnr, wni, h, n, npml, yr, yi, xr, xi, n_threads);
    
    if (xi_alloc) {  mxFree(xi); }

    return;
}
