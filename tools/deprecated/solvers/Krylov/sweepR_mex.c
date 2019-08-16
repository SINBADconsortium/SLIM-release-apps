/*
 * CARP/Kaczmarz sweep on diagonally banded matrix. The essential loop is:
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
 *		dir 		- if dir > 0, go through matrix rows in ascending order.
 * 			  		  if dir < 0, go through matrix rows in descending order.
 * 		n_threads	- OPTIONAL argument to control the number of execution
 * 		threads solving CARP blocks in parallel. The number of threads can also
 * 		be defined via an environment variable (OMP_NUM_THREADS), but this
 * 		optional argument takes precedence. The default number of threads is
 * 		one. Take care if using more than one MATLAB worker per node: each
 * 		MATLAB worker will use OMP_NUM_THREADS, so if there are four workers on
 * 		a node, there will be 4 x OMP_NUM_THREADS parallel CARP sweeps.
 *
 * Author: Art Petrenko, Tristan van Leeuwen
 *         Seismic Laboratory for Imaging and Modeling
 *         Department of Earth, Ocean, and Atmosperic Sciences
 *         The University of British Columbia
 *         
 * Date: July, 2014
 
 * You may use this code only under the conditions and terms of the
 * license contained in the file LICENSE provided with this source
 * code. If you do not agree to these terms you may not use this
 * software.
*/

#include <stdlib.h> /* for getenv */
#include <stddef.h> /* for size_t type */
#include <string.h> /* for memcpy */
#include <pthread.h> /* for threading */

/* The following section allows this file to compile on Mac OS X 10.8.5. Pass
 * the flag -DARCH_MACI64 to the compiler to activate it. */
#ifdef ARCH_MACI64
#include "pthread_barrier.h"
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
	int dir;
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
	long start_row, end_row, haloWidth, ncol, ny, nx;
	/* Rr and Ri are pointers to short fat ncol-by-N matrices */
	double *Rr, *Ri, *yr, *yi, *br, *bi;
	long *idx;
	double w;
	int dir;
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
	haloWidth	= thread_args->haloWidth;
	main_diagonal_offset = thread_args->main_diagonal_offset;
	Rr 			= thread_args->Rr;
	Ri 			= thread_args->Ri;
	idx 		= thread_args->idx;
	yr 			= thread_args->yr;
	yi 			= thread_args->yi;
	br 			= thread_args->br;
	bi 			= thread_args->bi;
	w 			= thread_args->w;
	dir 		= thread_args->dir;

	offset = (start_row == 0 ? 0 : haloWidth - main_diagonal_offset);

	/* Kaczmarz sweep on one row block */
	for(long i = (dir > 0 ? start_row : end_row-1);
		dir > 0 ? i<end_row : i>=start_row;
		dir > 0 ? i++ : i--)
	{

		if (0 <= i + main_diagonal_offset && i + main_diagonal_offset < nx){
			cr = br[i + main_diagonal_offset];
			ci = bi[i + main_diagonal_offset];
		}
		else{
		  //error(1,0,"Discovery of whether the iterate vector is haloed failed.");
		}
		/* First loop over non-zero row elements calculates inner product
		 * of matrix row and CARP iterate */
		for(long j=0, k; j<ncol;j++){
			/* i + idx[j] is the column index for the full Helmholtz matrix.
			 * k is the index into the vector representing the CARP iterate
			 * of the given block. */
			k = i + idx[j] - start_row + offset;
			if(0<=k && k<ny){
				cr -= Rr[i*ncol + j]*yr[k] - Ri[i*ncol + j]*yi[k];
				ci -= Rr[i*ncol + j]*yi[k] + Ri[i*ncol + j]*yr[k];
			}
		}
		/* Second loop over non-zero row elements updates CARP iterate */
		cr *= w;
		ci *= w;
		for(long j=0, k; j<ncol;j++){
			k = i + idx[j] - start_row + offset;
			if(0<=k && k<ny){
				yr[k] +=   cr*Rr[i*ncol + j] + ci*Ri[i*ncol + j];
				yi[k] +=  -cr*Ri[i*ncol + j] + ci*Rr[i*ncol + j];
			}
		}
	}

	return NULL;
}

void *average_halos(void *thread_args_void)
{	
	struct average_data_t *thread_args;
	double *copy_src_real, *copy_dst_real;
	double *copy_src_imag, *copy_dst_imag;
	double *halo_1_real, *halo_2_real;
	double *halo_1_imag, *halo_2_imag;
	double *halo_dst_real, *halo_dst_imag;
	size_t n_to_copy, n_in_halo;

	/* Assign local pointer to data locations in shared memory */
	thread_args 	= (struct average_data_t *) thread_args_void;
	copy_src_real	= thread_args->copy_src_real;
	copy_dst_real	= thread_args->copy_dst_real;
	copy_src_imag	= thread_args->copy_src_imag;
	copy_dst_imag	= thread_args->copy_dst_imag;
	n_to_copy		= thread_args->n_to_copy;
	halo_1_real     = thread_args->halo_1_real;
	halo_2_real		= thread_args->halo_2_real;
	halo_dst_real	= thread_args->halo_dst_real;
	halo_1_imag		= thread_args->halo_1_imag;
	halo_2_imag		= thread_args->halo_2_imag;
	halo_dst_imag   = thread_args->halo_dst_imag;
	n_in_halo		= thread_args->n_in_halo;

	/* Copy the non-halo parts of the domain block directly to output array */
	memcpy((void *)copy_dst_real, (void *)copy_src_real, sizeof(double)*n_to_copy);
	memcpy((void *)copy_dst_imag, (void *)copy_src_imag, sizeof(double)*n_to_copy);

	/* Average the halo parts of the domain block and copy to output array.
	 * NOTE: this assumes domain blocks overlap only with their nearest
	 * neighbour segments. */
	for (long i = 0; i < n_in_halo; i++){
		halo_dst_real[i] = (halo_1_real[i] + halo_2_real[i]) / 2;
		halo_dst_imag[i] = (halo_1_imag[i] + halo_2_imag[i]) / 2;
	}

	return NULL;
}

void *sweep_and_average(void *thread_data_void)
{
	struct thread_data_t *thread_data = (struct thread_data_t *) thread_data_void;

	/* Copy the initial guess into the working arrays */
	memcpy((void *)thread_data->copy_init_guess_data.copy_dst_real, (void *)thread_data->copy_init_guess_data.copy_src_real, sizeof(double)*thread_data->copy_init_guess_data.n_to_copy);
	memcpy((void *)thread_data->copy_init_guess_data.copy_dst_imag, (void *)thread_data->copy_init_guess_data.copy_src_imag, sizeof(double)*thread_data->copy_init_guess_data.n_to_copy);
	do_sweep((void *) &(thread_data->sweep_data));

	pthread_barrier_wait(thread_data->barrier);

	average_halos((void *) &(thread_data->average_data));

	return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* structs to hold all arguments to each thread in one variable */
	struct sweep_data_t **thread_args_sweep = NULL;
	struct average_data_t **thread_args_average = NULL;
	struct thread_data_t **thread_data = NULL;
	struct copy_init_guess_data_t **copy_init_guess_data = NULL;
	pthread_barrier_t barrier_after_sweep_before_halo_average;

	mwSize ncol, nx;
	ptrdiff_t ncolBlas, idxIncBlas = 1, maxIdxLoc = 0;
	mwSize haloWidth;
	mwSize n_threads = 1;
	char *n_threads_str = NULL;
	double *Rr,*Ri,*idxd = NULL,*xr,*xi,*br,*bi,*yr,*yi;
	long *idx = NULL;
    double w = 0;
	int dir = 1;
	mwSize N=1, numGridPointsPerBlock, main_diagonal_offset;

	/* Flags that are set if memory is allocated within the MEX file */
	int Ri_alloc=0, xi_alloc=0, bi_alloc=0;

	/* a return code flag and segment demarcation arrays */ 
	mwSize i_thread;
	mwSize *seg_bounds_hi, *seg_bounds_mid, *seg_bounds_row, *seg_bounds_lo;

	/* Threading variables */
	pthread_t *threadIDs = NULL;
	pthread_attr_t attr;

	/* Arrays to hold (overlapping) segments of yr and yi, a pair for each
	 * thread */
	double **yr_seg = NULL, **yi_seg = NULL;

	/* Allow worker threads to join back to main thread once they are done */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	/* Read input arguments; initialize complex part to zero if input is real. */
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
        xi = mxCalloc(nx,sizeof(double));
		xi_alloc = 1;
    }
	br  = mxGetPr(prhs[3]);
    if(mxIsComplex(prhs[3])){
        bi  = mxGetPi(prhs[3]);
    }
    else{
        bi = mxCalloc(nx,sizeof(double));
		bi_alloc = 1;
    }
	if (mxGetM(prhs[3]) != nx){
	    mexPrintf('%d %d \n',mxGetM(prhs[3]),nx);
		mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:NumElements", 
				"The number of elements in the iterate and right hand side vectors must be equal.");
	}
	w   = mxGetScalar(prhs[4]);
	dir = lrint(mxGetScalar(prhs[5]));

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

	/* Allocate the final output vector */
	plhs[0]  = mxCreateDoubleMatrix(nx, 1, mxCOMPLEX); 
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
	for (mwSize i=0; i < ncol; i++){
		idx[i] = lrint(idxd[i]); 
	}

	/* Compute (half) halo width. Remember that BLAS routines like idamax
	 * return FORTRAN-style indices. */
	maxIdxLoc = idamax_(&ncolBlas, idxd, &idxIncBlas) - 1;
	haloWidth = (mwSize)labs(idx[maxIdxLoc]); 

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
	numGridPointsPerBlock = N / n_threads;
	seg_bounds_hi  = (mwSize *)mxCalloc(n_threads+1,sizeof(mwSize));
	seg_bounds_mid = (mwSize *)mxCalloc(n_threads+1,sizeof(mwSize));
	seg_bounds_lo  = (mwSize *)mxCalloc(n_threads+1,sizeof(mwSize));
	seg_bounds_row = (mwSize *)mxCalloc(n_threads+1,sizeof(mwSize));
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
		for (mwSize i = 1; i <= ncol/2; i++){
			if (idx[ncol/2 + i] - main_diagonal_offset != -(idx[ncol/2 - i] - main_diagonal_offset)){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:DiagsNotSymmetric",
						"Input iterate vector appears to be haloed but the pattern of non-zero diagonals in the system matrix is not symmetric.");
			}
		}
	}
	for (i_thread=0; i_thread<n_threads; i_thread++){
		/* First domain block */
		if (i_thread==0){
			seg_bounds_hi[i_thread]  = 0;
			seg_bounds_mid[i_thread] = 0;
			seg_bounds_row[i_thread] = 0;
			seg_bounds_lo[i_thread]  = 0;
		}
		/* Other domain blocks */
		else {
			if (i_thread <= N % n_threads) {
				seg_bounds_mid[i_thread] = (numGridPointsPerBlock+1)*i_thread + main_diagonal_offset;
				seg_bounds_row[i_thread] = (numGridPointsPerBlock+1)*i_thread;
				seg_bounds_hi[i_thread]  = (numGridPointsPerBlock+1)*i_thread + haloWidth + main_diagonal_offset;	
				seg_bounds_lo[i_thread]  = (numGridPointsPerBlock+1)*i_thread - haloWidth + main_diagonal_offset;
			}
			else{
				seg_bounds_mid[i_thread] = (N % n_threads) + numGridPointsPerBlock*i_thread + main_diagonal_offset;
				seg_bounds_row[i_thread] = (N % n_threads) + numGridPointsPerBlock*i_thread;
				seg_bounds_hi[i_thread]  = (N % n_threads) + numGridPointsPerBlock*i_thread + haloWidth + main_diagonal_offset;	
				seg_bounds_lo[i_thread]  = (N % n_threads) + numGridPointsPerBlock*i_thread - haloWidth + main_diagonal_offset;
			}
			/* Check that halos do not overlap each other */
			if (seg_bounds_lo[i_thread] < seg_bounds_hi[i_thread-1]){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:TooManyThreads",
						"Too many threads; non-adjacent domain blocks share nearest-neighbour grid points.");
			}
		}
	}
	seg_bounds_lo[n_threads]  = nx;
	seg_bounds_hi[n_threads]  = nx;
	seg_bounds_mid[n_threads] = nx;
	seg_bounds_row[n_threads] = N;

	/* Allocate pointers to segments, to thread arguments and thread IDs */
	thread_args_sweep = (struct sweep_data_t **)mxCalloc(n_threads, sizeof(struct sweep_data_t *));
	if (n_threads > 1) {
		/* Set up a barrier for synchronization of all threads save the master that
		 * executes mexFunction. Note that strictly speaking, only threads working
		 * on domain blocks that share halos need to synchronize with each other.
		 * */
		if (pthread_barrier_init(&barrier_after_sweep_before_halo_average, NULL, n_threads)){
			mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:pthreads",
					"Could not initialize pthread barrier.");
		}

		thread_data = (struct thread_data_t **)mxCalloc(n_threads, sizeof(struct thread_data_t *));
		thread_args_average = (struct average_data_t **)mxCalloc(n_threads, sizeof(struct average_data_t *));
		copy_init_guess_data = (struct copy_init_guess_data_t **)mxCalloc(n_threads, sizeof(struct copy_init_guess_data_t *));
		yr_seg = (double **)mxCalloc(n_threads,sizeof(double *)); 
		yi_seg = (double **)mxCalloc(n_threads,sizeof(double *));
		threadIDs = (pthread_t *)mxCalloc(n_threads,sizeof(pthread_t));
	}
	for (i_thread=0; i_thread<n_threads; i_thread++){
		/* Allocate the segments for each thread */
		if (n_threads > 1) {
			yr_seg[i_thread] = (double *)mxCalloc((seg_bounds_hi[i_thread+1]-seg_bounds_lo[i_thread]),sizeof(double));
			yi_seg[i_thread] = (double *)mxCalloc((seg_bounds_hi[i_thread+1]-seg_bounds_lo[i_thread]),sizeof(double));
			/* Check that segments were allocated correctly */
			if (yr_seg[i_thread]==NULL || yi_seg[i_thread]==NULL){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:pthreadsOutOfMemory",
						"Could not allocate memory for thread computational variables.");
			}
			copy_init_guess_data[i_thread] = (struct copy_init_guess_data_t *)mxCalloc(1,sizeof(struct copy_init_guess_data_t));
			copy_init_guess_data[i_thread]->copy_src_real = &xr[seg_bounds_lo[i_thread]];
			copy_init_guess_data[i_thread]->copy_src_imag = &xi[seg_bounds_lo[i_thread]];
			copy_init_guess_data[i_thread]->copy_dst_real = yr_seg[i_thread];
			copy_init_guess_data[i_thread]->copy_dst_imag = yi_seg[i_thread];
			copy_init_guess_data[i_thread]->n_to_copy = (size_t) (seg_bounds_hi[i_thread+1] - seg_bounds_lo[i_thread]);
		}

		/* Set thread arguments */
		thread_args_sweep[i_thread] = (struct sweep_data_t *)mxCalloc(1,sizeof(struct sweep_data_t));
		thread_args_sweep[i_thread]->start_row 	= seg_bounds_row[i_thread];
		thread_args_sweep[i_thread]->end_row   	= seg_bounds_row[i_thread+1];
		thread_args_sweep[i_thread]->ncol 	 	= ncol;
		thread_args_sweep[i_thread]->ny 	 	= seg_bounds_hi[i_thread+1]-seg_bounds_lo[i_thread];
		thread_args_sweep[i_thread]->nx			= nx;
		thread_args_sweep[i_thread]->haloWidth 	= haloWidth*(i_thread ? 1 : 0);
		thread_args_sweep[i_thread]->main_diagonal_offset = main_diagonal_offset;
		thread_args_sweep[i_thread]->Rr  		= Rr;
		thread_args_sweep[i_thread]->Ri  		= Ri;
		thread_args_sweep[i_thread]->idx 		= idx;
		if (n_threads > 1){
			thread_args_sweep[i_thread]->yr 	= yr_seg[i_thread];
			thread_args_sweep[i_thread]->yi 	= yi_seg[i_thread];
		}
		else{
			thread_args_sweep[i_thread]->yr	 	= yr;
			thread_args_sweep[i_thread]->yi	 	= yi;
		}
		thread_args_sweep[i_thread]->br 		= br;
		thread_args_sweep[i_thread]->bi 		= bi;
		thread_args_sweep[i_thread]->w 		 	= w;
		thread_args_sweep[i_thread]->dir 		= dir;

		if (n_threads > 1){
			/* Set the arguments for the averaging threads. Note that each
			 * thread is responsible for averaging the low end of its address
			 * range. The middle of its address range is copied, while the high
			 * end of its address range is left for the next thread. */
			thread_args_average[i_thread] = (struct average_data_t *)mxCalloc(1,sizeof(struct average_data_t));
			thread_args_average[i_thread]->copy_src_real = &(yr_seg[i_thread][seg_bounds_hi[i_thread] - seg_bounds_lo[i_thread]]);
			thread_args_average[i_thread]->copy_dst_real = &(yr[seg_bounds_hi[i_thread]]);
			thread_args_average[i_thread]->copy_src_imag = &(yi_seg[i_thread][seg_bounds_hi[i_thread] - seg_bounds_lo[i_thread]]);
			thread_args_average[i_thread]->copy_dst_imag = &(yi[seg_bounds_hi[i_thread]]);
			thread_args_average[i_thread]->n_to_copy 	 = (size_t) (seg_bounds_lo[i_thread+1] - seg_bounds_hi[i_thread]);
			thread_args_average[i_thread]->halo_1_real   = &(yr_seg[i_thread-1 >= 0 ? i_thread-1 : 0][seg_bounds_lo[i_thread] - seg_bounds_lo[i_thread-1 >= 0 ? i_thread-1 : 0]]);
			thread_args_average[i_thread]->halo_2_real	 = yr_seg[i_thread];
			thread_args_average[i_thread]->halo_dst_real = &(yr[seg_bounds_lo[i_thread]]);
			thread_args_average[i_thread]->halo_1_imag   = &(yi_seg[i_thread-1 >= 0 ? i_thread-1 : 0][seg_bounds_lo[i_thread] - seg_bounds_lo[i_thread-1 >= 0 ? i_thread-1 : 0]]);
			thread_args_average[i_thread]->halo_2_imag	 = yi_seg[i_thread];
			thread_args_average[i_thread]->halo_dst_imag = &(yi[seg_bounds_lo[i_thread]]);
			thread_args_average[i_thread]->n_in_halo 	 = (size_t) (seg_bounds_hi[i_thread] - seg_bounds_lo[i_thread]);

			thread_data[i_thread] = (struct thread_data_t *)mxCalloc(1,sizeof(struct thread_data_t));
			thread_data[i_thread]->copy_init_guess_data = *copy_init_guess_data[i_thread];
			thread_data[i_thread]->sweep_data = *thread_args_sweep[i_thread];
			thread_data[i_thread]->average_data = *thread_args_average[i_thread];
			thread_data[i_thread]->barrier = &barrier_after_sweep_before_halo_average;
		}
	}

	if (n_threads == 1) {
		/* Set the initial guess directly in the output array too */
		memcpy((void *)yr, (void *)xr, sizeof(double)*nx);
		memcpy((void *)yi, (void *)xi, sizeof(double)*nx);
		do_sweep((void *)thread_args_sweep[0]);
	}
	else{
		/* Sweep and average, in separate worker threads */
		for (i_thread=0; i_thread<n_threads; i_thread++){
			if (pthread_create(&threadIDs[i_thread], &attr, sweep_and_average, (void *)thread_data[i_thread])){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:pthreadCreate", 
						"pthread_create returned non-zero.");
			}
		}

		/* Wait for threads to finish */
		for (i_thread=0; i_thread<n_threads; i_thread++){
			if (pthread_join(threadIDs[i_thread], NULL)){
				mexErrMsgIdAndTxt("SLIM_release_apps:tools:algorithms:ThreeDFreqModeling:sweepR_mex:pthreadCreate", 
						"pthread_join returned non-zero.");
			}
		}
	}

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
	mxFree(seg_bounds_hi);
	mxFree(seg_bounds_mid);
	mxFree(seg_bounds_row);
	mxFree(seg_bounds_lo);
	if (n_threads > 1) {
		mxFree(threadIDs);
	}
	for (i_thread=0; i_thread<n_threads; i_thread++){
		mxFree(thread_args_sweep[i_thread]);
		if (n_threads > 1){
			mxFree(thread_args_average[i_thread]);
			mxFree(thread_data[i_thread]);
			mxFree(copy_init_guess_data[i_thread]);
			mxFree(yr_seg[i_thread]);
			mxFree(yi_seg[i_thread]);
		}
	}
	if (n_threads > 1){
		pthread_barrier_destroy(&barrier_after_sweep_before_halo_average);
		mxFree(thread_args_average);
		mxFree(thread_data);
		mxFree(copy_init_guess_data);
		mxFree(yr_seg);
		mxFree(yi_seg);
	}
	mxFree(idx);
	mxFree(thread_args_sweep);
	pthread_attr_destroy(&attr);

	/* Don't think I need pthread_exit() here, because pthread_join is called above */
	return;
}
