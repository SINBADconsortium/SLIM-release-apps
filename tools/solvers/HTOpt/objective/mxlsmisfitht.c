#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <stdlib.h>
#include "string.h"
#include <stdio.h>

#define exp2(a) (1 << (a))
#define mod2(a,d) ((a) % exp2(d))
#define tree2lin(h,k) (exp2(h-1) - 1 + k)
#define idx1d(a,b,dim1) (a + (b)*(dim1))
#define idx1d3(a,b,c,dim1,dim2) (a + (b)*(dim1) + (c)*(dim2)*(dim1))
#define Bidx1d(h,k) (exp2(h)+k-1)
#define max(a,b) (a > b ? a : b)

/*
  Implementation of Algorithm 2 in 'Optimization on the Hierarchical Tucker manifold - applications to tensor completion' - C. Da Silva and F. J. Herrmann.
 */

void
printRow(const double* A)
{
    printf("%3.2f %3.2f %3.2f %3.2f %3.2f\n",A[0],A[1],A[2],A[3],A[4]);
}
void
printU(const double* U)
{
    for(int i=0; i<5; i++)
    {
        for(int j=0; j<10; j++)
        {
            printf("%3.3f ",U[i +j*5]);
        }
        printf("\n");
    }
}
void
printTens(const double* B,int m,int n, int p)
{
  for(int k=0; k<p; k++)
  {
    printf("k: %d \n",k);
    for(int i=0; i<m; i++)
    {
      for(int j=0; j<n; j++)
      {
	printf("%3.4f ",B[i + j*m + k*m*n]);
      }
      printf("\n");
    }
    printf("\n\n");
  }
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define U prhs[0]
#define B prhs[1]
#define I prhs[2]
#define b prhs[3]
#define Uidx1d prhs[4]
#define stopidx prhs[5]

#define dU plhs[1]
#define dB plhs[2]
    
    unsigned int * I_p;
    int numSamples, d, depth;
    mwSize * dims;
    double ** U_p, *** B_p;
    mxArray * Utmp, * dUtmp;
    double * Utmp_p, * b_p;
    unsigned int * Uidx1d_p;
    double ** dU_p, *** dB_p;
    size_t stop_idx,K,nnodes;
    
    int Iidx, Uidx1d_pos, U1d_pos,endj;
    int this_idx, left_idx, right_idx;
    int kl, kr, kt;
    size_t * ranks;
    double f = 0;
    double Uval, r, M;
    int i,j,ii,jj,kk,i_samp;
    
    numSamples = mxGetN(I);
    d = mxGetM(I);
    
    depth = (int)ceil(log2(d)) + 1;
    stop_idx = (int) *mxGetPr(stopidx);
    //M = sqrt(numSamples);
    M = 1;
    nnodes = exp2(depth) + 2*( d % exp2(depth-1) ) - 1;
    
    //Get relevant pointers/arrays from inputs
    I_p = (unsigned int *) mxGetData(I);
    Uidx1d_p = (unsigned int *)mxGetData(Uidx1d);
    b_p = (double * )mxGetPr(b);
    U_p = (double**)mxCalloc(d,sizeof(double *));
    B_p = (double***)mxCalloc(depth-1,sizeof(double**));
    dims = (mwSize*)mxCalloc(d,sizeof(mwSize));
    ranks = (mwSize * )mxCalloc(nnodes,sizeof(mwSize));
    const mwSize * tmpSize;
    ranks[0] = 1;
    K = 1;
    for(i=1;i<mxGetM(B);i++)
    {
      tmpSize = mxGetDimensions(mxGetCell(B,i));
      ranks[i] = tmpSize[2];
      K = max(K,ranks[i]);
    }

    for(i=0; i<d; i++)
    {
        U_p[i] = mxGetPr(mxGetCell(U,i));
        dims[i] = mxGetN(mxGetCell(U,i));
	ranks[Uidx1d_p[i]] = mxGetM(mxGetCell(U,i));
	K = max(K,ranks[Uidx1d_p[i]]);
    }
    for(i=0; i<depth-1; i++)
    {
        if (i==depth-2)
        {
            B_p[i] = (double**) mxCalloc(stop_idx,sizeof(double*));
            for(j=0; j<stop_idx; j++)
            {
                B_p[i][j] = mxGetPr(mxGetCell(B,Bidx1d(i,j)));
            }
        }
        else
        {
            B_p[i] = (double**) mxCalloc(exp2(i),sizeof(double*));
            for(j=0; j<exp2(i); j++)
            {
                B_p[i][j] = mxGetPr(mxGetCell(B,Bidx1d(i,j)));
            }
        }
    }
    
    //Temporary array for the HT tree
    Utmp = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetM(Utmp,K); mxSetN(Utmp,nnodes);
    mxSetData(Utmp, mxMalloc(sizeof(double) * K * nnodes) );
    Utmp_p = mxGetPr(Utmp);
    
    if(nlhs == 1)
    {                
        for(i_samp=0; i_samp<numSamples;i_samp++)
        {
            //Copy values of U{i} in to temporary array
            for(i=0; i<d; i++)
            {
                Iidx = I_p[ idx1d(i,i_samp,d) ] - 1;
                //U{i}(I(i_samp),1:k)
		kt = ranks[Uidx1d_p[i]];
                U1d_pos = idx1d(0,Iidx,kt);
                
                //Position in the 1D tree
                Uidx1d_pos = idx1d(0,Uidx1d_p[i],K);
                
                //Copy U{i}(:,I(i_samp)) in to the temporary array
                memcpy(&( Utmp_p[ Uidx1d_pos ] ), &( U_p[i][ U1d_pos ] ), sizeof(double)* kt);
            }
         
	    //Compute U
	    for( i=depth-2; i >= 0; i--)
	    {
	      if (i==depth-2)
		endj = stop_idx;
	      else
		endj = exp2(i);
	      for( j=0; j < endj; j++)
	      {
		this_idx = exp2(i) + j-1;
		left_idx = exp2(i+1) + 2*(j)-1;
		right_idx = left_idx+1;
		kl = ranks[left_idx];
		kr = ranks[right_idx];
		kt = ranks[this_idx];
		if (i==0)
		{
		  Uval = 0;
		  for(ii=0; ii<kl; ii++)
		  {
		    for(jj=0; jj<kr; jj++)
		    {
		      Uval += Utmp_p[ idx1d(ii,left_idx,K) ] * B_p[0][0][ idx1d(ii,jj,kl) ] * Utmp_p[ idx1d(jj,right_idx,K) ];
		    }
		  }
		}
		else
		{
		  for(kk=0; kk<kt; kk++)
		  {
		    Uval = 0;
		    for(ii=0; ii<kr; ii++)
		    {
		      for(jj=0; jj<kl; jj++)
		      {
			Uval += Utmp_p[ idx1d(jj,left_idx,K) ] * Utmp_p[ idx1d(ii,right_idx,K) ] * B_p[i][j][ idx1d3(jj,ii,kk,kl,kr) ];
		      }
		    }                                
		    Utmp_p[ idx1d(kk,this_idx,K) ] = Uval;
		  }
		}
	      }
	    }
	    //Residual
	    r = (Uval - b_p[i_samp])/M;
	    f += 0.5*r*r;
            
            //Assign output variables
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            double * f_out = mxGetPr(plhs[0]);
            f_out[0] = f;
        }                
    }
    else if (nlhs == 3)
    {
      //Temporary array to hold tree information
      dUtmp = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetM(dUtmp,K); mxSetN(dUtmp,nnodes);
      mxSetData(dUtmp, mxMalloc(sizeof(double) * K * nnodes) );
      double * dUtmp_p = mxGetPr(dUtmp);
      double dULval, dURval, dBval, s, t, u, v,w;
      dB_p = (double***)mxCalloc(depth-1,sizeof(double**));
      dU_p = (double**)mxCalloc(d,sizeof(double*));
      for(i=0; i<d; i++)
      {
	kt = ranks[ Uidx1d_p[i] ];
	dU_p[i] = (double*) mxCalloc(dims[i]*kt,sizeof(double));
      }
	      
      for(i=0; i<depth-1; i++)
      {
	if (i==depth-2)                    
	  endj = stop_idx;                    
	else                    
	  endj = exp2(i);
	
	dB_p[i] = (double**) mxCalloc(endj,sizeof(double**));
	for(j=0; j<endj; j++)
	{
	  this_idx = exp2(i) + j-1;
	  left_idx = exp2(i+1) + 2*(j)-1;
	  right_idx = left_idx+1;
	  kl = ranks[left_idx];
	  kr = ranks[right_idx];
	  kt = ranks[this_idx];
	  dB_p[i][j] = (double*)mxCalloc(kl*kr*kt,sizeof(double));
	}	 
      }
      for(i_samp=0; i_samp<numSamples;i_samp++)
      {
	//Copy values of U{i} in to temporary array
	for(i=0; i<d; i++)
	{
	  Iidx = I_p[ idx1d(i,i_samp,d) ] - 1;
	  //U{i}(I(i_samp),1:k)
	  kt = ranks[ Uidx1d_p[i] ];

	  U1d_pos = idx1d(0,Iidx,kt);
                
	  //Position in the 1D tree
	  Uidx1d_pos = idx1d(0,Uidx1d_p[i],K);
                
	  //Copy U{i}(:,I(i_samp)) in to the temporary array
	  memcpy(&( Utmp_p[ Uidx1d_pos ] ), &( U_p[i][ U1d_pos ] ), sizeof(double)*kt );
	}
                          
	//Compute U
	for( i=depth-2; i >= 0; i--)
	{
	  if (i==depth-2)                    
	    endj = stop_idx;                    
	  else                    
	    endj = exp2(i);
                    
	  for( j=0; j < endj; j++)
	  {
	    this_idx = exp2(i) + j-1;
	    left_idx = exp2(i+1) + 2*(j)-1;
	    right_idx = left_idx+1;
	    kl = ranks[left_idx];
	    kr = ranks[right_idx];
	    kt = ranks[this_idx];
	    if (i==0)
	    {
	      Uval = 0;
	      for(ii=0; ii<kl; ii++)
	      {
		for(jj=0; jj<kr; jj++)
		{
		  Uval += Utmp_p[ idx1d(ii,left_idx,K) ] * B_p[0][0][ idx1d(ii,jj,kl) ] * Utmp_p[ idx1d(jj,right_idx,K) ];
		}
	      }
	    }
	    else
	    {
	      for(kk=0; kk<kt; kk++)
	      {
		Uval = 0;
		for(ii=0; ii<kr; ii++)
		{
		  for(jj=0; jj<kl; jj++)
		  {
		    Uval += Utmp_p[ idx1d(jj,left_idx,K) ] * Utmp_p[ idx1d(ii,right_idx,K) ] * B_p[i][j][ idx1d3(jj,ii,kk,kl,kr) ];
		  }
		}                                
		Utmp_p[ idx1d(kk,this_idx,K) ] = Uval;
	      }
	    }
	  }
	}
	//Residual
	r = (Uval - b_p[i_samp])/M;
	f += 0.5*r*r;
                
	for( i=0; i<depth-1; i++)
	{
	  if (i==depth-2)                    
	    endj = stop_idx;                    
	  else                    
	    endj = exp2(i);
	  for( j=0; j < endj; j++ )
	  {
	    this_idx = exp2(i) + j-1;
	    left_idx = exp2(i+1) + 2*(j)-1;
	    right_idx = left_idx+1;
	    kl = ranks[left_idx];
	    kr = ranks[right_idx];
	    kt = ranks[this_idx];
	    if (i==0)
	    {
	      for(ii=0; ii < kr; ii++)
	      {
		dURval = 0;
		for(jj=0; jj < kl; jj++)
		{                                    
		  dB_p[0][0][ idx1d(jj,ii,kl) ] += r * Utmp_p[ idx1d(jj,left_idx,K) ] * Utmp_p[ idx1d(ii,right_idx,K) ];
		  dURval += Utmp_p[ idx1d(jj,left_idx,K) ] * B_p[0][0][ idx1d(jj,ii,kl) ];		                    
		}	
		dUtmp_p[ idx1d(ii,right_idx,K) ] = r*dURval;
	      }
	      for(ii=0; ii<kl; ii++)
	      {
		dULval = 0; 
		for(jj=0; jj<kr; jj++)
		{
		  dULval += Utmp_p[ idx1d(jj,right_idx,K) ] * B_p[0][0][ idx1d(ii,jj,kl) ];
		}
		dUtmp_p[ idx1d(ii,left_idx,K) ] = r*dULval;	
	      }
	    }
	    else
	    {
	      //ii indexes over kt dimension
	      for(ii=0; ii < kt; ii++)
	      {	
		s = dUtmp_p[ idx1d(ii,this_idx,K) ];
		for(jj=0; jj< kr; jj++)
		{                                                                       		  
		  for(kk=0; kk<kl; kk++)
		  {
		    dB_p[i][j][ idx1d3(kk,jj,ii,kl,kr) ] += Utmp_p[ idx1d(kk,left_idx,K) ] * s * Utmp_p[ idx1d(jj,right_idx,K) ];
		  }
		}		
	      }  
	      //ii indexes over kr dimension
	      for(ii=0; ii<kr; ii++)
	      {
		dURval = 0;  
		for(jj=0; jj<kt; jj++)
		{
		  w = dUtmp_p[ idx1d(jj,this_idx,K) ];
		  for(kk=0; kk<kl; kk++)
		  {
		    dURval += w * Utmp_p[ idx1d(kk,left_idx,K) ] * B_p[i][j][ idx1d3(kk,ii,jj,kl,kr) ];
		  }
		}
		dUtmp_p[ idx1d(ii,right_idx,K) ] = dURval;
	      }
	      //ii indexes over kl dimension
	      for(ii=0; ii<kl; ii++)
	      {
		dULval = 0;
		for(jj=0; jj<kt; jj++)
		{
		  w = dUtmp_p[ idx1d(jj,this_idx,K) ];
		  for(kk=0; kk<kr; kk++)
		  {
		    dULval += w * Utmp_p[ idx1d(kk,right_idx,K) ] * B_p[i][j][ idx1d3(ii,kk,jj,kl,kr) ];
		  }
		}
		dUtmp_p[ idx1d(ii,left_idx,K) ] = dULval;
	      }	     
	    }
	  }                    
	}
	//Accumulate variables at output
	for(i=0; i<d; i++)
	{
	  Iidx = I_p[ idx1d(i,i_samp,d) ] - 1;
	  kt = ranks[ Uidx1d_p[i] ];
	  for(j=0;j<kt; j++)
	  {
	    dU_p[i][ idx1d(j,Iidx,kt) ] += dUtmp_p[ idx1d(j,Uidx1d_p[i],K) ];
	  }	  	 	    
	} 
      }
      //Assign output variables
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      double * f_out = mxGetPr(plhs[0]);
      f_out[0] = f;
      
      mxArray * tmpOut;
      mwSize emptysize[3];
      emptysize[0] = 0; emptysize[1] = 0; emptysize[2] = 0;
      mwSize bsize[3];
      bsize[0] = K; bsize[1] = K; bsize[2] = K;
      mwIndex ij[2];

      dU = mxCreateCellMatrix(1,d);
      for(i=0; i<d; i++)
      {
	tmpOut = mxCreateDoubleMatrix(0,0,mxREAL);
	kt = ranks[ Uidx1d_p[i] ];
	mxSetM(tmpOut,kt);
	mxSetN(tmpOut,dims[i]);
	mxSetData(tmpOut,dU_p[i]);
	mxSetCell(dU,i,tmpOut);
      }
      
      dB = mxCreateCellMatrix(mxGetM(B),mxGetN(B));
      int numB = 0;
      for(i=0; i<depth-1;i++)
      {
	if (i==depth-2)                    
	  endj = stop_idx;                    
	else                    
	  endj = exp2(i);
	for( j=0; j < endj; j++ )
	{
	  this_idx = exp2(i) + j-1;
	  left_idx = exp2(i+1) + 2*(j)-1;
	  right_idx = left_idx+1;
	  kl = ranks[left_idx];
	  kr = ranks[right_idx];
	  kt = ranks[this_idx];
	  bsize[0] = kl; bsize[1] = kr; bsize[2] = kt;
	  tmpOut = mxCreateNumericArray(3,emptysize,mxDOUBLE_CLASS,0);
	  mxSetDimensions(tmpOut,bsize,3);
	  mxSetData(tmpOut,dB_p[i][j]);
	  ij[0] = i; ij[1] = j;
	  mxSetCell(dB,numB,tmpOut);
	  numB++;
	}		
      }
      
      //Clear derivative specific memory
      //dU_p[i] and dB_p[i][j] aren't freed because they're the results returned to matlab
      //These are just the pointer variables allocated in this program
      for(i=0; i<depth-1; i++)
      {
	mxFree(dB_p[i]);       
      }
      mxFree(dB_p);
      mxFree(dU_p);
    }     
    
    //Free created memory
    mxDestroyArray(Utmp);
    mxFree(U_p);
    mxFree(dims);
    for(i=0; i<depth-1; i++)
    {
        mxFree(B_p[i]);
    }
    mxFree(B_p);
    mxFree(ranks);
    
    
}
