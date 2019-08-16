#include <stdlib.h> /* for getenv */
#include <stddef.h> /* for size_t type */
#include <string.h> /* for memcpy */
#include <omp.h>

#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "Helm3d_27pt.h"

#ifdef DERIV
#undef DERIV
#endif

/*
 *
 * 27pt stencil Helmholtz matrix-vector product.  
 * 
 * 
 * Usage:
 *  y = Helm3dkaczswp_mex(wn,h,n,nmpl,omega,x,b);
 *
 * Input:
 *   wn       - wavenumber evaluated on pml grid
 *   h        - [hx,hy,hz] spacings in each direction
 *   n        - [nx,ny,nz] points in each dimension, including pml
 *   npml     - 2 x 3 matrix denoting the number of pml points in each dimension                    
 *   omega    - relaxation parameter, in (0,2)
 *   x        - nx*ny*nz x 1 vector
 *   b        - right hand side, nx*ny*nz x 1 vector
 *
 * Output:
 *   y      - Kaczmarz sweep applied to x
 *  
 * Curt Da Silva, 2015
 *
 * To compile, run:
 *   mex -O -largeArrayDims Helm3dkaczswp_mex.c -DDEFINEUNIX -lmwblas CFLAGS="\$CFLAGS -fopenmp -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 */

#define WN prhs[0]
#define H prhs[1]
#define N prhs[2]
#define NPML prhs[3]
#define OMEGA prhs[4]
#define X prhs[5]
#define B prhs[6]
#define OUTPUT plhs[0]

void init_zero(double *x, int n){
    int i;
    
    for(i=0;i<n;i++){
        x[i] = 0;
    }
}

void do_sweep( double * wnr, double * wni, double * h, double * n, double * npml, 
	       double *yr, double *yi, double *xr, double *xi, double * br, double * bi, 
	       const double omega, int P) {      	 
    int i,j,k,kout,t,s;

    
    int nx = (int)n[0]; int ny = (int)n[1]; int nz = (int)n[2];
    
    int npmlx_lo = (int)npml[0]; int npmlx_hi = (int)npml[1]; 
    int npmly_lo = (int)npml[2]; int npmly_hi = (int)npml[3]; 
    int npmlz_lo = (int)npml[4]; int npmlz_hi = (int)npml[5];

    coef_consts c = compute_coef_consts(h);    

    int pmlz_alloc = 0; int pmly_alloc = 0; int pmlx_alloc = 0;

    double complex coef[27];     
    double complex x[27]; 
    double complex d;
    double row_norm_sq;
    
    wn_type wn_window[27];

    pml_info p;    
    pml_adj_info padj;
	
    double complex atx;            
    p.x_hasL = 0; p.x_hasR = 1;
    p.y_hasL = 0; p.y_hasR = 1;
    p.z_hasL = 0; p.z_hasR = 1;
#pragma omp parallel for schedule(static) private(i,j,k,kout,t,wn_window,coef,atx,row_norm_sq,d,x) firstprivate(p,padj,pmlz_alloc,pmly_alloc,pmlx_alloc) num_threads(P)
    for(s=0; s<P; s++)
    {
	xyzloop_updown( 
	    // Cache a window of the wavenumber around the current point 
	    load_wn_nbrhood(wn_window,wnr,wni,i,j,k,nx,ny,nz,p);
	    
	    // Get coefficients 
	    get_coefs(coef,wn_window, c, p, padj);		
	    
	    // Cache a window of the wavefield around the current point 
	    load_nbrhoodc(x,yr,yi,i,j,k,nx,ny,nz,s,p);
	    
	    kout = IDX1D4(i,j,k,s,nx,ny,nz);
	    atx = 0.0 + 0.0*I;
	    row_norm_sq = 0;
	    for(t=0; t<27; t++){
		atx += coef[t] * x[t];
		row_norm_sq += creal(conj(coef[t])*coef[t]);
	    }
	    
	    d = CMPLX(br[kout],bi[kout]);
	    d = omega*(d - atx)/row_norm_sq;
	    
	    for(t=0; t<27; t++){
		coef[t] = d*conj(coef[t]);
	    }
	    nbrhood_update(coef,yr,yi,i,j,k,nx,ny,nz,s,p);	
	    )
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    mwSize n_threads = 1;
        
    double * wnr, *wni, *h, *n, *npml, *yr, *yi, *br, *bi, *xr, *xi;     
    double omega;
    char *n_threads_str = NULL;       
    int bi_alloc = 0;
    int xi_alloc = 0;
    int P;
    wnr   = mxGetPr(WN);
    if(mxIsComplex(WN))
    {
	wni   = mxGetPi(WN);
    }
    else
    {
	wni = NULL;
    }
                            
    if(mxGetM(N) != 3 && mxGetN(N) != 3)
    {
	mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:Helm3dkaczswp_mex:Nsize",
			  "n must be a 3-length vector");
    }

    h     = mxGetPr(H);
    n     = mxGetPr(N);
    omega = mxGetScalar(OMEGA);
    npml  = mxGetPr(NPML);   

    int numel = (int)n[0]*n[1]*n[2];

    if(mxGetM(X) != numel) {
        mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:Helm3dkaczswp_mex:Xsize",
			  "x must have nx*ny*nz elements");
    }
    P = mxGetN(X);
    if(mxGetM(B) != numel) {
        mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:Helm3dkaczswp_mex:Bsize",
			  "b must have nx*ny*nz elements");
    }

    xr  = mxGetPr(X);
    if(mxIsComplex(X)){
        xi  = mxGetPi(X);
    }
    else{
        xi = mxCalloc(numel*P,sizeof(double));
	xi_alloc = 1;
    }
    
    
    br = mxGetPr(B);
    if(mxIsComplex(B)){ bi = mxGetPi(B); }
    else
    {
	bi = mxCalloc(numel*P,sizeof(double));
	bi_alloc = 1;
    }    
    

     /* define output vector y and initialize with input vector x.*/
    OUTPUT   = mxCreateDoubleMatrix(numel, P, mxCOMPLEX);    
    yr       = mxGetPr(OUTPUT);
    yi       = mxGetPi(OUTPUT);
    
    memcpy( yr, xr, sizeof(double)*numel*P );
    memcpy( yi, xi, sizeof(double)*numel*P );

    do_sweep( wnr, wni, h, n, npml, yr, yi, xr, xi, br, bi,omega,P );
    
    if (bi_alloc) {  mxFree(bi); }
    if (xi_alloc) {  mxFree(xi); }
    
    return;
}
