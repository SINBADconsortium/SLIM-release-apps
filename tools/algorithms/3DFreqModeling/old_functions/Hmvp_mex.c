#include "math.h"
#include "mex.h"   
#include "matrix.h"
#include <assert.h>

#define COEF prhs[0]
#define IDX prhs[1]
#define X prhs[2]
#define NTHREADS prhs[3]
#define OUTPUT plhs[0]

/*
IMPORTANT: This function is a modified version of sweep_mex. The comments are
still the original from the sweeps though. It means that most comments will not
make any sense. Someone should fix this - I have already fixed way too many 
things!! - Lago

H multiplication for 3D penalty method
for each row i
	x = x + w*(b(i) - A(i,:)*x)*A(i,:)'
end
The matrix is given in banded storage format, where each column of the array S
stores a band of the matrix such that A(i,i+idx(j)) = S(i,j).

use:
	y = sweepR_mex(S.',idx,x,b,w,dir,n)
NOTE: this function expects S.', not S.
*/

void init_zero(double *x, int n){
   int i;

   for(i=0;i<n;i++){
      x[i] = 0;
   }
   }

   void init_array(double *y, double *x, int n){
   int i;

   for(i=0;i<n;i++){
      y[i] = x[i];
   }
}

/* Sr and Si are pointers to a short fat ncol x nrow matrix holding the bands
 * of the Helmholtz operator as rows. */
void do_Hmvp(int nrow, int ncol, int nx, double *Sr, double *Si, double *idx, double *yr, double *yi, double *xr, double *xi){
   int i,j,k,l,m,n;
   /*actual sweep, forward direction*/
   /*i loops over the columns */
   l  = (int)idx[0];
   l  = abs(l);
   for(i=0;i<nrow;i++){
      /* first calculate the inner products*/
      /* Note that idx[13] refers to the main diagonal of the Helmholtz
         * matrix. The index into the idx array will change with the
         * stencil used. */
      /*j loops over bands (rows) */
      m = i + l;
      n = i*ncol;
      for(j=0;j<ncol;j++){
         /*k is the column index for the full Helmholtz matrix A */
         k = m + (int)idx[j];
         yr[i] += Sr[n + j]*xr[k] - Si[n + j]*xi[k];
         yi[i] += Sr[n + j]*xi[k] + Si[n + j]*xr[k];
      }
   }
/* same, in reverse direction*/
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
   double *Sr,*Si,*xr,*xi,*br,*bi,*yr,*yi;
   double *idx, w;
   int dir, nrow, ncol, nx;
   int Si_alloc=0, xi_alloc=0, bi_alloc=0;
   int i;

   /* read input, initialize complex part to zero if input is real.*/
   nrow = mxGetN(COEF);
   ncol = mxGetM(COEF);
   Sr  = mxGetPr(COEF);
   if(mxIsComplex(COEF)){
      Si  = mxGetPi(COEF);
   }
   else{
      Si = mxCalloc(nrow*ncol,sizeof(double));
      init_zero(Si,nrow*ncol);
      Si_alloc = 1;
   }
   idx = mxGetPr(IDX);	
   nx  = mxGetM(X);
   xr  = mxGetPr(X);
   if(mxIsComplex(X)){
      xi  = mxGetPi(X);
   }
   else{
   xi = mxCalloc(nx,sizeof(double));
      init_zero(xi,nx);
      xi_alloc = 1;
   }

   /* define output vector y and initialize with input vector x.*/
   OUTPUT   = mxCreateDoubleMatrix(nrow, 1, mxCOMPLEX); 
   yr       = mxGetPr(OUTPUT);
   yi       = mxGetPi(OUTPUT);
   init_zero(yr,nrow);
   init_zero(yi,nrow);

   /*sweep*/
   do_Hmvp(nrow, ncol, nx, Sr,Si,idx,yr,yi,xr,xi);

   if (Si_alloc){
      mxFree(Si);
   }
   if (xi_alloc){
      mxFree(xi);
   }

   return;
}
