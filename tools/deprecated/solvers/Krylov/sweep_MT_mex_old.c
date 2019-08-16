/*
 * Kaczmarz down-up sweep on diagonally banded matrix. The essential loop is:
 *
 * for each row i
 *     x = x + w*(b(i) - A(i,:)*x)*A(i,:)'
 * end
 *
 * The matrix is given in band storage format, where each row (stored
 * contiguously in memory) of the array R stores a diagonal of the matrix with
 * offset idx(j), such that 
 *
 * 	   A(i,i+idx(j)) = R(i,j)
 *
 * use (from MATLAB):
 * 		y = sweepR_mex(R,idx,x,b,w,dir)
 *
 * 		R			- matrix of diagonals of matrix A
 *		idx			- offsets of diagonals
 *		x			- initial guess
 * 		b 			- right hand side (source)
 *		w 			- relaxation parameter (0 <= w <= 2)
 *		ns 		    - number of RHS
 * 		n_threads	- OPTIONAL argument to control the number of execution
 * 		threads solving CARP blocks in parallel. The number of threads can also
 * 		be defined via an environment variable (OMP_NUM_THREADS), but this
 * 		optional argument takes precedence. The default number of threads is
 * 		one. Take care if using more than one MATLAB worker per node: each
 * 		MATLAB worker will use OMP_NUM_THREADS, so if there are four workers on
 * 		a node, there will be 4 x OMP_NUM_THREADS parallel CARP sweeps.
 *  compile :
 *	   mex -largeArrayDims sweep_MT_mex.c -DDEFINEUNIX -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 * Author: Mathias Louboutin from Art Petrenko sweepR_mex.c
 *         Seismic Laboratory for Imaging and Modeling
 *         Department of Earth, Ocean, and Atmosperic Sciences
 *         The University of British Columbia
 *         
 * Date: March, 2015
 
 * You may use this code only under the conditions and terms of the
 * license contained in the file LICENSE provided with this source
 * code. If you do not agree to these terms you may not use this
 * software.
*/

#include <stdlib.h> /* for getenv */
#include <stddef.h> /* for size_t type */
#include <string.h> /* for memcpy */
#include <pthread.h> /* for threading */
#include <omp.h>

/* The following section allows this file to compile on Mac OS X 10.8.5. Pass
 * the flag -DARCH_MACI64 to the compiler to activate it. */
#ifdef ARCH_MACI64
#include <mach/error.h>
typedef wchar_t char16_t;
#else /* not ARCH_MACI64 */
#include <error.h>
#endif /* ARCH_MACI64 */

#include <math.h>
#include <mex.h>   
#include <matrix.h>
#include <blas.h>

struct copy_init_guess_data_t 
{
	double *copy_src_real, *copy_dst_real;
	double *copy_src_imag, *copy_dst_imag;
	long n_to_copy;
};

struct sweep_data_t 
{
	long start_row, end_row, ncol, ny, nx, haloWidth, main_diagonal_offset;
	double *Rr, *Ri, *yr, *yi, *br, *bi;
	long *idx;
	double w;
	int dir, n_threads, ns;
};

struct average_data_t 
{
	double *copy_src_real, *copy_dst_real;
	double *copy_src_imag, *copy_dst_imag;
	double *halo_1_real, *halo_2_real;
	double *halo_1_imag, *halo_2_imag;
	double *halo_dst_real, *halo_dst_imag;
	long n_to_copy, n_in_halo;
};

struct thread_data_t 
{
	struct copy_init_guess_data_t copy_init_guess_data;
	struct sweep_data_t sweep_data;
	struct average_data_t average_data;
	pthread_barrier_t *barrier;
};

void *do_sweep(void *thread_args_void)
{
    struct sweep_data_t *thread_args;
	/* Variables contained in thread_args_void struct */
	long start_row, end_row, ncol, ny, nx;
	/* Rr and Ri are pointers to short fat ncol-by-N matrices */
	double *Rr, *Ri, *yr, *yi, *br, *bi;
	long *idx;
	double w;
	int n_threads,ns;
	/* Temporary storage variables */
	double cr = 0, ci = 0;
	long offset, main_diagonal_offset;
	
	/* Assign local pointers to data locations in shared memory */
	thread_args = (struct sweep_data_t *) thread_args_void;
	start_row 	= thread_args->start_row;
	end_row 	= thread_args->end_row;
	ncol 		= thread_args->ncol;
	ny 			= thread_args->ny;
	nx			= thread_args->nx;
	main_diagonal_offset = thread_args->main_diagonal_offset;
	Rr 			= thread_args->Rr;
	Ri 			= thread_args->Ri;
	idx 		= thread_args->idx;
	yr 			= thread_args->yr;
	yi 			= thread_args->yi;
	br 			= thread_args->br;
	bi 			= thread_args->bi;
	w 			= thread_args->w;
	n_threads 	= thread_args->n_threads;
	ns 			= thread_args->ns;
	offset = (start_row == 0 ? 0 : - main_diagonal_offset);
	
	long i;
	int s;
	long j;
	long k;
	long indj;
	long toto;

	#pragma omp parallel for schedule(static,1) private(k,indj,i,j,cr,ci) num_threads(n_threads)
	for(s=0 ; s<ns; s++)
	{
		/* Kaczmarz sweep on one row block */
		for(i = start_row ; i<end_row ;i++ )
		{

			if (0 <= i + main_diagonal_offset && i + main_diagonal_offset < nx){
				cr = br[i + main_diagonal_offset+s*nx];
				ci = bi[i + main_diagonal_offset+nx*s];
			}
			else{
				error(1,0,"Discovery of whether the iterate vector is haloed failed.");
			}
			/* First loop over non-zero row elements calculates inner product
			 * of matrix row and CARP iterate */
			long diff = i - start_row + offset;
			long icol=i*ncol;

			for(j=0 ; j<ncol;j++)			
			{
				/* i + idx[j] is the column index for the full Helmholtz matrix.
				 * k is the index into the vector representing the CARP iterate
				 * of the given block. */
				k = diff + idx[j];
				indj=icol + j;
				if(0<=k && k<ny)
				{
					cr -= Rr[indj]*yr[k+s*nx] - Ri[indj]*yi[k+s*nx];
					ci -= Rr[indj]*yi[k+s*nx] + Ri[indj]*yr[k+s*nx];

				}
			}

			cr*=w;
			ci*=w;
			/* Second loop over non-zero row elements updates Karkzmarz iterate */

			for(j=0 ; j<ncol;j++)
			{
				k = diff + idx[j];
				indj=icol + j;
				if(0<=k && k<ny)
				{
					yr[k+s*nx] +=   cr*Rr[indj] + ci*Ri[indj];
					yi[k+s*nx] +=  -cr*Ri[indj] + ci*Rr[indj];

				}
			}


		}
		cr = 0;
		ci = 0;
		for(i = end_row-1 ; i>start_row-1 ;i-- )
		{

			if (0 <= i + main_diagonal_offset && i + main_diagonal_offset < nx){
				cr = br[i + main_diagonal_offset+s*nx];
				ci = bi[i + main_diagonal_offset+nx*s];
			}
			else{
				error(1,0,"Discovery of whether the iterate vector is haloed failed.");
			}
			/* First loop over non-zero row elements calculates inner product
			 * of matrix row and CARP iterate */
			long diff = i - start_row + offset;
			long icol=i*ncol;

			for(j=0 ; j<ncol;j++)			
			{
				/* i + idx[j] is the column index for the full Helmholtz matrix.
				 * k is the index into the vector representing the CARP iterate
				 * of the given block. */
				k = diff + idx[j];
				indj=icol + j;
				if(0<=k && k<ny)
				{
					cr -= Rr[indj]*yr[k+s*nx] - Ri[indj]*yi[k+s*nx];
					ci -= Rr[indj]*yi[k+s*nx] + Ri[indj]*yr[k+s*nx];

				}
			}

			cr*=w;
			ci*=w;
			/* Second loop over non-zero row elements updates Karkzmarz iterate */

			for(j=0 ; j<ncol;j++)
			{
				k = diff + idx[j];
				indj=icol + j;
				if(0<=k && k<ny)
				{
					yr[k+s*nx] +=   cr*Rr[indj] + ci*Ri[indj];
					yi[k+s*nx] +=  -cr*Ri[indj] + ci*Rr[indj];

				}
			}


		}
	}
	
	

	return NULL;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* structs to hold all arguments to each thread in one variable */
	struct sweep_data_t **thread_args_sweep = NULL;
	struct copy_init_guess_data_t **copy_init_guess_data = NULL;
	
	mwSize ncol, nx;
	ptrdiff_t ncolBlas, idxIncBlas = 1, maxIdxLoc = 0;

	double *Rr,*Ri,*idxd = NULL,*xr,*xi,*br,*bi,*yr,*yi;
	long *idx = NULL;
    double w = 0;
	int ns;
	mwSize n_threads = 1;
	char *n_threads_str = NULL;
	mwSize N=1, numGridPointsPerBlock, main_diagonal_offset;

	/* Flags that are set if memory is allocated within the MEX file */
	int Ri_alloc=0, xi_alloc=0, bi_alloc=0;
	mwSize *seg_bounds_hi, *seg_bounds_mid, *seg_bounds_row, *seg_bounds_lo;
	
	/* Read input arguments; initialize complex part to zero if input is real. */
	ns = lrint(mxGetScalar(prhs[5]));
	N = mxGetN(prhs[0]);
    ncol = mxGetM(prhs[0]);
	ncolBlas = (ptrdiff_t)ncol;
	Rr  = mxGetPr(prhs[0]);
    if(mxIsComplex(prhs[0])){
        Ri  = mxGetPi(prhs[0]);
    }
    else{
        Ri = mxCalloc(N*ncol,sizeof(double));
		Ri_alloc = 1;
    }
	idxd = mxGetPr(prhs[1]);	
    nx  = mxGetM(prhs[2]);
	xr  = mxGetPr(prhs[2]);
    if(mxIsComplex(prhs[2])){
        xi  = mxGetPi(prhs[2]);
    }
    else{
        xi = mxCalloc(nx*ns,sizeof(double));
		xi_alloc = 1;
    }
	br  = mxGetPr(prhs[3]);
    if(mxIsComplex(prhs[3])){
        bi  = mxGetPi(prhs[3]);
    }
    else{
        bi = mxCalloc(nx*ns,sizeof(double));
		bi_alloc = 1;
    }
    if (mxGetM(prhs[3]) != nx){
    mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:NumElements", 
				"The number of elements in the iterate and right hand side vectors must be equal.");
	}



	/* Allocate the final output vector */
	plhs[0]  = mxCreateDoubleMatrix(nx,ns, mxCOMPLEX); 
	yr       = mxGetPr(plhs[0]);
    yi       = mxGetPi(plhs[0]);

	/* Check to make sure memory was allocated correctly */
	if (Rr==NULL || Ri==NULL || idxd==NULL || xr==NULL || xi==NULL || br==NULL || bi==NULL ||
	    yr==NULL || yi==NULL){
		mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:OutOfMemory", 
				"Could not allocate memory for main computational variables.");
	}

	if ((idx = (long *) mxCalloc(ncol, sizeof(long))) == NULL){
		mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:OutOfMemory",
				"Could not allocate memory for main computational variables.");
	}
	mwSize i;
	for (i=0; i < ncol; i++){
		idx[i] = lrint(idxd[i]); 
	}
	
	/* The default value for the number of threads can be overridden by an
	 * environment variable. */
	n_threads_str = getenv("OMP_NUM_THREADS");
	if (n_threads_str == NULL){
		n_threads = 1;
	}
	else{
		n_threads = strtol(n_threads_str, NULL, 10);
		if(n_threads < 1){
			n_threads = 1;
		}
	}
	/* The environment variable can in turn be overridden by an optional
	 * argument to the mexFunction. */ 
	if (nrhs >= 7){
		if(1 <= lrint(mxGetScalar(prhs[6]))){
			n_threads = lrint(mxGetScalar(prhs[6]));
		}
	}

	/* printf("Using %d threads \n",n_threads);
	/* Partition the iterate vector into blocks. Note that the below
	 * partitioning scheme is slighlty different from that in pCARPCG.m in this
	 * directory. The partitioning scheme of pCARPCG corresponds to
	 * distributing a three dimensional array with dimensions given by n
	 * according to Matlab's codistributor1d.defaultPartition(n(3)), and then
	 * vectorizing it. The partition scheme of the present file instead uses
	 * Matlab's codistributor1d.defaultPartion(prod(n)). In other words,
	 * pCARPCG divides the iterate into blocks along the slow dimension,
	 * whereas sweepR_mex.c does not take dimensionality into account, only the
	 * total number of gridpoints. This is done to avoid needing an extra input
	 * parameter with the the system dimensions. The seg_bounds_hi, _lo and
	 * _mid arrays contain indices into non-haloed vectors, while the
	 * seg_bounds_row array contains indices to the rows of the system matrix.
	 *
		 * yr_seg[i_thread-1]        overlap      yr_seg[i_thread]
	 * ------------------------|-----|-----|-------------------------------
	 *         .----------------^     |     ^-------------------.
	 *    seg_bounds_lo[i_thread], seg_bounds_mid[i_thread], seg_bounds_hi[i_thread]
	 */
	numGridPointsPerBlock = N;
	seg_bounds_hi  = (mwSize *)mxCalloc(2,sizeof(mwSize));
	seg_bounds_mid = (mwSize *)mxCalloc(2,sizeof(mwSize));
	seg_bounds_lo  = (mwSize *)mxCalloc(2,sizeof(mwSize));
	seg_bounds_row = (mwSize *)mxCalloc(2,sizeof(mwSize));
	if (N == nx){
		main_diagonal_offset = 0;
	}
	else{
		/* The vector is haloed. We are only able to correctly process matrices
		 * with a non-zero main diagonal and symmetric off-main diagonal offsets. */
		if (ncol % 2 != 1){
			mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:EvenNumberOfDiags",
					"Input iterate vector appears to be haloed but there is an even number of non-zero diagonals in the system matrix.");
		}
		main_diagonal_offset = idx[ncol/2];
		mwSize i;
		for (i = 1; i <= ncol/2; i++){
			if (idx[ncol/2 + i] - main_diagonal_offset != -(idx[ncol/2 - i] - main_diagonal_offset)){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:DiagsNotSymmetric",
						"Input iterate vector appears to be haloed but the pattern of non-zero diagonals in the system matrix is not symmetric.");
			}
		}
	}
	
	seg_bounds_hi[0]  = 0;
	seg_bounds_mid[0] = 0;
	seg_bounds_row[0] = 0;
	seg_bounds_lo[0]  = 0;
	seg_bounds_lo[1]  = nx;
	seg_bounds_hi[1]  = nx;
	seg_bounds_mid[1] = nx;
	seg_bounds_row[1] = N;
	
	thread_args_sweep = (struct sweep_data_t **)mxCalloc(1, sizeof(struct sweep_data_t *));

	/* Set thread arguments */
	thread_args_sweep[0] = (struct sweep_data_t *)mxCalloc(1,sizeof(struct sweep_data_t));
	thread_args_sweep[0]->start_row 	= seg_bounds_row[0];
	thread_args_sweep[0]->end_row   	= seg_bounds_row[1];
	thread_args_sweep[0]->ncol 	 	= ncol;
	thread_args_sweep[0]->ny 	 	= seg_bounds_hi[1]-seg_bounds_lo[0];
	thread_args_sweep[0]->nx			= nx;
	thread_args_sweep[0]->main_diagonal_offset = main_diagonal_offset;
	thread_args_sweep[0]->Rr  		= Rr;
	thread_args_sweep[0]->Ri  		= Ri;
	thread_args_sweep[0]->idx 		= idx;
	thread_args_sweep[0]->yr	 	= yr;
	thread_args_sweep[0]->yi	 	= yi;
	thread_args_sweep[0]->br 		= br;
	thread_args_sweep[0]->bi 		= bi;
	thread_args_sweep[0]->w 		 = w;
	thread_args_sweep[0]->n_threads = n_threads;
	thread_args_sweep[0]->ns = ns;
	/* Set the initial guess directly in the output array too */
	memcpy((void *)yr, (void *)xr, sizeof(double)*nx*ns);
	memcpy((void *)yi, (void *)xi, sizeof(double)*nx*ns);
	do_sweep((void *)thread_args_sweep[0]);

	/* Free memory if it was allocated within the MEX file. */
	if (Ri_alloc){
		mxFree(Ri);
	}
	if (xi_alloc){
		mxFree(xi);
	}
	if (bi_alloc){
		mxFree(bi);
	}

	mxFree(idx);
	mxFree(thread_args_sweep);

	/* Don't think I need pthread_exit() here, because pthread_join is called above */
	return;
}
