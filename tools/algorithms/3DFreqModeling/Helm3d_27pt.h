#include "math.h"
#include <complex.h>

#ifdef SINGLE_PRECISION
#define numeric float;
#else
#define numeric double;
#endif

typedef double _Complex cmplx_type;

#ifdef WNCMPLX
typedef cmplx_type wn_type;
#else
typedef double wn_type;
#endif

#define PI M_PI

#define CMPLX(x, y) ((double complex)((double)(x) + _Complex_I * (double)(y)))
#define IDX1D(a,b,m,n) ((a) + (b)*(m))
#define IDX1D3(a,b,c,m,n,p) ((a) + (b)*(m) + (c)*(m)*(n))
#define IDX1D4(a,b,c,d,m,n,p) ((a) + (b)*(m) + (c)*(m)*(n) + (d)*(m)*(n)*(p))

/*
void print_c(const double complex c)
{
    mexPrintf("%3.15e + i%3.15e ",creal(c),cimag(c));
}

void print_nbrhdc(const double complex * c)
{
    for(int k=0; k<3; k++)
    {
	for(int j=0; j<3; j++)
	{
	    for(int i=0; i<3; i++)
	    {
		print_c(c[IDX1D3(i,j,k,3,3,3)]); mexPrintf("\n");
	    }
	    mexPrintf("\n");
	}
	mexPrintf("\n\n");
    }
}
void print_nbrhd(const double * c)
{
    for(int k=0; k<3; k++)
    {
	for(int j=0; j<3; j++)
	{
	    for(int i=0; i<3; i++)
	    {
		mexPrintf("%3.15e  ",c[IDX1D3(i,j,k,3,3,3)]); mexPrintf("\n");
	    }
	    mexPrintf("\n");
	}
	mexPrintf("\n\n");
    }
}
*/

#define xyzloop_updown(...) for(k=0; k<nz; k++)				    \
    {									    \
	pmlzfunc(&p,&padj,k,nz,npmlz_lo,npmlz_hi,pmlz_alloc);		    \
	pmlz_alloc = 1;							    \
	p.z_hasL = k>0; p.z_hasR = k < nz-1;				    \
	for(j=0; j<ny; j++)						    \
	{								    \
	    p.y_hasL = j > 0; p.y_hasR = j < ny-1;			    \
	    pmlyfunc(&p,&padj,j,ny,npmly_lo,npmly_hi,pmly_alloc);	    \
	    pmly_alloc = 1;						    \
	    for(i=0; i<nx; i++)						    \
	    {								    \
		p.x_hasL = i > 0; p.x_hasR = i < nx-1;			    \
		pmlxfunc(&p,&padj,i,nx,npmlx_lo,npmlx_hi,pmlx_alloc,true);  \
		pmlx_alloc = 1;						    \
		__VA_ARGS__						    \
	    }							            \
	    pmlx_alloc = 0;						    \
	}								    \
	pmly_alloc = 0;	     						    \
    }									    \
    pmlz_alloc = 0;  							    \
    for(k=nz-1; k>=0; k--)						    \
    {								            \
        pmlzfunc(&p,&padj,k,nz,npmlz_lo,npmlz_hi,pmlz_alloc);               \
	pmlz_alloc = 1;				                            \
	p.z_hasL = k>0; p.z_hasR = k < nz-1;                                \
        for(j=ny-1; j>=0; j--)	                                            \
	{								    \
	    p.y_hasL = j > 0; p.y_hasR = j < ny-1;			    \
	    pmlyfunc(&p,&padj,j,ny,npmly_lo,npmly_hi,pmly_alloc);	    \
	    pmly_alloc = 1;						    \
	    for(i=nx-1; i>=0; i--)					    \
	    {			                                            \
		p.x_hasL = i > 0; p.x_hasR = i < nx-1;			    \
		pmlxfunc(&p,&padj,i,nx,npmlx_lo,npmlx_hi,pmlx_alloc,false); \
		pmlx_alloc = 1;						    \
		__VA_ARGS__						    \
	    }							            \
	    pmlx_alloc = 0;						    \
	}								    \
        pmly_alloc = 0;							    \
    }									    \
    pmlz_alloc = 0;

#define xyzloop_up(...) for(k=0; k<nz; k++)				    \
    {									    \
	pmlzfunc(&p,&padj,k,nz,npmlz_lo,npmlz_hi,pmlz_alloc);		    \
	pmlz_alloc = 1;							    \
	p.z_hasL = k > 0; p.z_hasR = k < nz-1;				    \
	for(j=0; j<ny; j++)						    \
	{								    \
	    p.y_hasL = j > 0; p.y_hasR = j < ny-1;			    \
	    pmlyfunc(&p,&padj,j,ny,npmly_lo,npmly_hi,pmly_alloc);	    \
	    pmly_alloc = 1;						    \
	    for(i=0; i<nx; i++)						    \
	    {								    \
		p.x_hasL = i > 0; p.x_hasR = i < nx-1;			    \
		pmlxfunc(&p,&padj,i,nx,npmlx_lo,npmlx_hi,pmlx_alloc,1);  \
		pmlx_alloc = 1;						    \
		__VA_ARGS__						    \
	    }							            \
	    pmlx_alloc = 0;						    \
	}								    \
	pmly_alloc = 0;	     						    \
    }

#define xyzloop_down(...)  for(k=nz-1; k>=0; k--)  			    \
    {									    \
	pmlzfunc(&p,&padj,k,nz,npmlz_lo,npmlz_hi,pmlz_alloc);		    \
	pmlz_alloc = 1;							    \
        p.z_hasL = k > 0; p.z_hasR = k < nz-1;				    \
	for(j=ny-1; j>=0; j--)					            \
	{								    \
	    pmlyfunc(&p,&padj,j,ny,npmly_lo,npmly_hi,pmly_alloc);	    \
	    p.y_hasL = j > 0; p.y_hasR = j < ny-1;			    \
	    pmly_alloc = 1;						    \
	    for(i=nx-1; i>=0; i--)         				    \
	    {								    \
	        p.x_hasL = i > 0; p.x_hasR = i < nx-1;			    \
		pmlxfunc(&p,&padj,i,nx,npmlx_lo,npmlx_hi,pmlx_alloc,0); \
		pmlx_alloc = 1;						    \
		__VA_ARGS__						    \
	    }							            \
	    pmlx_alloc = 0;						    \
	}								    \
	pmly_alloc = 0;	     						    \
    }

/*
  The current coordinate in each direction is denoted as N
  +1/-1 in the corresponding coordinate is denoted as P/M, respectively
  So, e.g., 
  NNN - current point (x,y,z)
  MNN - (x-1,y,z)
  PNN - (x+1,y,z)
  MNP - (x-1,y,z+1)
  PMP - (x+1,y-1,z+1)
  etc.
*/        
int const MMM = 0;
int const NMM = 1;
int const PMM = 2;
int const MNM = 3;
int const NNM = 4;
int const PNM = 5;
int const MPM = 6;
int const NPM = 7;
int const PPM = 8;
int const MMN = 9;
int const NMN = 10;
int const PMN = 11;
int const MNN = 12;
int const NNN = 13;
int const PNN = 14;
int const MPN = 15;
int const NPN = 16;
int const PPN = 17;
int const MMP = 18;
int const NMP = 19;
int const PMP = 20;
int const MNP = 21;
int const NNP = 22;
int const PNP = 23;
int const MPP = 24;
int const NPP = 25;
int const PPP = 26;

/*
  Absolute constants related to the stencil
 */
double const W1 = (1.8395262e-5);
double const W2 = (0.296692333333333);
double const W3 = (0.027476150000000);
double const WM1 = (0.49649658);
double const WM2 = (0.075168750000000);
double const WM3 = (0.004373916666667);
double const WM4 = (5.690375e-07);

#ifdef DERIV
double const non_deriv_mode = 0;
#else
double const non_deriv_mode = 1;
#endif

/*
  Run time constants related to the stencil
 */
struct _coef_consts {
    const double wn_coef; 
    const double wn_xcoef;
    const double wn_ycoef;
    const double wn_zcoef;

    const double pmlx_coef;
    const double pmly_coef;
    const double pmlz_coef;

    const double xz_coef;
    const double xy_coef;
    const double yz_coef;

    const double W3A_2;    

};
typedef struct _coef_consts coef_consts;

/*
  Structs for holding information of pml function values in each direction
 */
struct _pml_info {
    double complex pmlx_lo; 
    double complex pmlx_hi;
    double complex pmlx;
    double complex pmly_lo;
    double complex pmly_hi;
    double complex pmly;
    double complex pmlz_lo;
    double complex pmlz_hi;
    double complex pmlz;
    int x_hasL;
    int x_hasR;
    int y_hasL;
    int y_hasR;
    int z_hasL;
    int z_hasR;
};
typedef struct _pml_info pml_info;

struct _pml_adj_info {
    double complex pmlx_lo_window[3];
    double complex pmlx_hi_window[3];
    double complex pmly_lo_window[3];
    double complex pmly_hi_window[3];
    double complex pmlz_lo_window[3];
    double complex pmlz_hi_window[3];
};
typedef struct _pml_adj_info pml_adj_info;


#define MMM_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasL && p.z_hasL,1) ? ( __VA_ARGS__ ) : 0
#define NMM_BDRY(...) __builtin_expect(p.y_hasL && p.z_hasL,1) ? ( __VA_ARGS__ ) : 0
#define PMM_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasL && p.z_hasL,1) ? ( __VA_ARGS__ ) : 0
#define MNM_BDRY(...) __builtin_expect(p.x_hasL && p.z_hasL,1) ? ( __VA_ARGS__ ) : 0
#define NNM_BDRY(...) __builtin_expect(p.z_hasL,1) ? ( __VA_ARGS__) : 0
#define PNM_BDRY(...) __builtin_expect(p.x_hasR && p.z_hasL,1) ? ( __VA_ARGS__) : 0

#define MPM_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasR && p.z_hasL,1) ? ( __VA_ARGS__) : 0
#define NPM_BDRY(...) __builtin_expect(p.y_hasR &&  p.z_hasL,1) ? ( __VA_ARGS__) : 0
#define PPM_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasR && p.z_hasL,1) ? ( __VA_ARGS__) : 0
#define MMN_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasL,1) ? ( __VA_ARGS__) : 0
#define NMN_BDRY(...) __builtin_expect(p.y_hasL,1) ? ( __VA_ARGS__) : 0
#define PMN_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasL,1) ? ( __VA_ARGS__) : 0
#define MNN_BDRY(...) __builtin_expect(p.x_hasL,1) ? ( __VA_ARGS__) : 0
#define PNN_BDRY(...) __builtin_expect(p.x_hasR,1) ? ( __VA_ARGS__) : 0
#define MPN_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasR,1) ? ( __VA_ARGS__) : 0
#define NPN_BDRY(...) __builtin_expect(p.y_hasR,1) ? ( __VA_ARGS__) : 0
#define PPN_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasR,1) ? ( __VA_ARGS__) : 0
#define MMP_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasL && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define NMP_BDRY(...) __builtin_expect(p.y_hasL && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define PMP_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasL &&  p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define MNP_BDRY(...) __builtin_expect(p.x_hasL && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define NNP_BDRY(...) __builtin_expect(p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define PNP_BDRY(...) __builtin_expect(p.x_hasR && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define MPP_BDRY(...) __builtin_expect(p.x_hasL && p.y_hasR && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define NPP_BDRY(...) __builtin_expect(p.y_hasR && p.z_hasR,1) ? ( __VA_ARGS__) : 0
#define PPP_BDRY(...) __builtin_expect(p.x_hasR && p.y_hasR && p.z_hasR,1) ? ( __VA_ARGS__) : 0

#define MMM_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasL && p.z_hasL,1)) { __VA_ARGS__ }
#define NMM_IF(...) if(__builtin_expect(p.y_hasL && p.z_hasL,1)){ __VA_ARGS__ }
#define PMM_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasL && p.z_hasL,1)){ __VA_ARGS__ }
#define MNM_IF(...) if(__builtin_expect(p.x_hasL && p.z_hasL,1)){ __VA_ARGS__ }
#define NNM_IF(...) if(__builtin_expect(p.z_hasL,1)) { __VA_ARGS__ }
#define PNM_IF(...) if(__builtin_expect(p.x_hasR && p.z_hasL,1)){ __VA_ARGS__ }

#define MPM_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasR && p.z_hasL,1)){ __VA_ARGS__ }
#define NPM_IF(...) if(__builtin_expect(p.y_hasR &&  p.z_hasL,1)) { __VA_ARGS__ }
#define PPM_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasR && p.z_hasL,1)) { __VA_ARGS__ }
#define MMN_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasL,1)){ __VA_ARGS__ }
#define NMN_IF(...) if(__builtin_expect(p.y_hasL,1)){ __VA_ARGS__ }
#define PMN_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasL,1)){ __VA_ARGS__ }
#define MNN_IF(...) if(__builtin_expect(p.x_hasL,1)){ __VA_ARGS__ }
#define PNN_IF(...) if(__builtin_expect(p.x_hasR,1)){ __VA_ARGS__ }
#define MPN_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasR,1)){ __VA_ARGS__ }
#define NPN_IF(...) if(__builtin_expect(p.y_hasR,1)){ __VA_ARGS__ }
#define PPN_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasR,1)){ __VA_ARGS__ }
#define MMP_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasL && p.z_hasR,1)){ __VA_ARGS__ }
#define NMP_IF(...) if(__builtin_expect(p.y_hasL && p.z_hasR,1)){ __VA_ARGS__ }
#define PMP_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasL &&  p.z_hasR,1)){ __VA_ARGS__ }
#define MNP_IF(...) if(__builtin_expect(p.x_hasL && p.z_hasR,1)){ __VA_ARGS__ }
#define NNP_IF(...) if(__builtin_expect(p.z_hasR,1)){ __VA_ARGS__ }
#define PNP_IF(...) if(__builtin_expect(p.x_hasR && p.z_hasR,1)){ __VA_ARGS__ }
#define MPP_IF(...) if(__builtin_expect(p.x_hasL && p.y_hasR && p.z_hasR,1)){ __VA_ARGS__ }
#define NPP_IF(...) if(__builtin_expect(p.y_hasR && p.z_hasR,1)){ __VA_ARGS__ }
#define PPP_IF(...) if(__builtin_expect(p.x_hasR && p.y_hasR && p.z_hasR,1)){ __VA_ARGS__ }

/******************************** 
   
     PML related functions

 ********************************/


inline double gamma_func_lower( double x, double nx, double npml_bot, double npml_top )
{
    return cos(PI*((x-1)* nx/(2*(nx+1)*npml_bot)));
}

inline double gamma_func_upper( double x, double nx, double npml_bot, double npml_top )
{
    return cos(PI*((1-(x-1)/(nx+1)) * nx/(2*npml_top)));
}

/*
 There are two functions associated to the pml, denoted pml_func_lower and pml_func_upper
 
 
 
 x - in the range [1,nx]
 nx - length of domain, including pml
 npml_bot - number of pml pts on the bottom of the domain
 npml_top - number of pml pts on the top of the domain
 */
inline double complex pml_func_lower( double x, double nx, double npml_bot, double npml_top )
{
    double gamma, gammap1;
    if (x <= npml_bot)
    {
	gammap1 = (x==npml_bot) ? 0 : gamma_func_lower(x+1,nx,npml_bot,npml_top);
	gamma = gamma_func_lower(x,nx,npml_bot,npml_top);	
    }
    else{
	if( (x > npml_bot) && (x <= nx+2-npml_bot) )
	{
	    gamma = 0; 
	    gammap1 = (x==nx+2-npml_bot) ? gamma_func_upper(x+1,nx,npml_bot,npml_top) : 0;
	}
	else
	{
	    gamma = gamma_func_upper(x,nx,npml_bot,npml_top);
	    gammap1 = gamma_func_upper(x+1,nx,npml_bot,npml_top);
	}
    }
    return (2.0 + 0.0*I)/( (2.0 - gammap1*(gammap1+gamma)) + (3.0*gammap1+gamma)*I);
}

inline double complex pml_func_upper( double x, double nx, double npml_bot, double npml_top )
{
    double gammap2, gammap1;
    if (x <= npml_bot)
    {
	gammap1 = (x==npml_bot) ? 0 : gamma_func_lower(x+1,nx,npml_bot,npml_top);
	gammap2 = (x >= npml_bot-1) ? 0 : gamma_func_lower(x+2,nx,npml_bot,npml_top);	
    }
    else{
	if( (x > npml_bot) && (x <= nx+2-npml_bot) )
	{
	    gammap2 = (x>=nx+1-npml_bot) ? gamma_func_upper(x+2,nx,npml_bot,npml_top) : 0; 
	    gammap1 = (x==nx+2-npml_bot) ? gamma_func_upper(x+1,nx,npml_bot,npml_top) : 0;
	}
	else
	{
	    gammap2 = gamma_func_upper(x+2,nx,npml_bot,npml_top);
	    gammap1 = gamma_func_upper(x+1,nx,npml_bot,npml_top);
	}
    }
    return (2.0 + 0.0*I)/( (2.0 - gammap1*(gammap1+gammap2)) + (3.0*gammap1+gammap2)*I);
}



inline void pmlzfunc(pml_info * p, pml_adj_info * padj, int k, int nz, int npmlz_lo, int npmlz_hi, int pmlz_alloc)
{
#ifdef ADJ
    padj->pmlz_lo_window[0] = pml_func_lower(k,nz,npmlz_lo,npmlz_hi);
    padj->pmlz_hi_window[0] = pml_func_upper(k,nz,npmlz_lo,npmlz_hi);
    padj->pmlz_lo_window[1] = pml_func_lower(k+1,nz,npmlz_lo,npmlz_hi);
    padj->pmlz_hi_window[1] = pml_func_upper(k+1,nz,npmlz_lo,npmlz_hi);
    padj->pmlz_lo_window[2] = pml_func_lower(k+2,nz,npmlz_lo,npmlz_hi);
    padj->pmlz_hi_window[2] = pml_func_upper(k+2,nz,npmlz_lo,npmlz_hi);
    
    p->pmlz = padj->pmlz_lo_window[1] + padj->pmlz_hi_window[1];
#else
    p->pmlz_lo = pml_func_lower(k+1,nz,npmlz_lo, npmlz_hi);
    p->pmlz_hi = pml_func_upper(k+1,nz,npmlz_lo, npmlz_hi);
    p->pmlz = p->pmlz_lo + p->pmlz_hi;
#endif
    
}

inline void pmlyfunc(pml_info * p, pml_adj_info * padj, int j, int ny, int npmly_lo, int npmly_hi, int pmly_alloc)
{
#ifdef ADJ    
    padj->pmly_lo_window[0] = pml_func_lower(j,ny,npmly_lo,npmly_hi);
    padj->pmly_hi_window[0] = pml_func_upper(j,ny,npmly_lo,npmly_hi);
    padj->pmly_lo_window[1] = pml_func_lower(j+1,ny,npmly_lo,npmly_hi);
    padj->pmly_hi_window[1] = pml_func_upper(j+1,ny,npmly_lo,npmly_hi);
    padj->pmly_lo_window[2] = pml_func_lower(j+2,ny,npmly_lo,npmly_hi);
    padj->pmly_hi_window[2] = pml_func_upper(j+2,ny,npmly_lo,npmly_hi);
    
    p->pmly = padj->pmly_lo_window[1] + padj->pmly_hi_window[1];
#else
    p->pmly_lo = pml_func_lower(j+1,ny,npmly_lo, npmly_hi);
    p->pmly_hi = pml_func_upper(j+1,ny,npmly_lo, npmly_hi);
    p->pmly = p->pmly_lo + p->pmly_hi;
#endif
    
}

inline void pmlxfunc(pml_info * p, pml_adj_info * padj, int i, int nx, int npmlx_lo, int npmlx_hi, int pmlx_alloc, int i_incr)
{
#ifdef ADJ       
    if(!pmlx_alloc)
    {
       padj->pmlx_lo_window[0] = pml_func_lower(i,nx,npmlx_lo,npmlx_hi);
       padj->pmlx_hi_window[0] = pml_func_upper(i,nx,npmlx_lo,npmlx_hi);
       padj->pmlx_lo_window[1] = pml_func_lower(i+1,nx,npmlx_lo,npmlx_hi);
       padj->pmlx_hi_window[1] = pml_func_upper(i+1,nx,npmlx_lo,npmlx_hi);
       padj->pmlx_lo_window[2] = pml_func_lower(i+2,nx,npmlx_lo,npmlx_hi);
       padj->pmlx_hi_window[2] = pml_func_upper(i+2,nx,npmlx_lo,npmlx_hi);
    }
    else
    {
	if(i_incr)
	{
	    padj->pmlx_lo_window[0] = padj->pmlx_lo_window[1]; 
	    padj->pmlx_lo_window[1] = padj->pmlx_lo_window[2]; 
	    padj->pmlx_lo_window[2] = pml_func_lower(i+2,nx,npmlx_lo,npmlx_hi);
	    padj->pmlx_hi_window[0] = padj->pmlx_hi_window[1]; 
	    padj->pmlx_hi_window[1] = padj->pmlx_hi_window[2]; 
	    padj->pmlx_hi_window[2] = pml_func_upper(i+2,nx,npmlx_lo,npmlx_hi);	    
	}
	else
	{
	    padj->pmlx_lo_window[2] = padj->pmlx_lo_window[1]; 
	    padj->pmlx_lo_window[1] = padj->pmlx_lo_window[0]; 
	    padj->pmlx_lo_window[0] = pml_func_lower(i,nx,npmlx_lo,npmlx_hi);
	    padj->pmlx_hi_window[2] = padj->pmlx_hi_window[1]; 
	    padj->pmlx_hi_window[1] = padj->pmlx_hi_window[0]; 
	    padj->pmlx_hi_window[0] = pml_func_upper(i,nx,npmlx_lo,npmlx_hi);
	}
    }
    p->pmlx = padj->pmlx_lo_window[1] + padj->pmlx_hi_window[1];
#else
    p->pmlx_lo = pml_func_lower(i+1,nx,npmlx_lo, npmlx_hi);
    p->pmlx_hi = pml_func_upper(i+1,nx,npmlx_lo, npmlx_hi);
    p->pmlx = p->pmlx_lo + p->pmlx_hi;
#endif    
}




/******************************** 
   
   Read/write 27pt neighbourhoods of real/complex values

 ********************************/


inline void load_nbrhoodc( double complex * x, const double * xr, const double * xi, int i, int j, int k, int nx, int ny, int nz, int s, pml_info p )
{    
    x[MMM] = MMM_BDRY(CMPLX( xr[ IDX1D4(i-1,j-1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j-1,k-1,s,nx,ny,nz) ] ));
    x[NMM] = NMM_BDRY(CMPLX( xr[ IDX1D4(i  ,j-1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j-1,k-1,s,nx,ny,nz) ] ));
    x[PMM] = PMM_BDRY(CMPLX( xr[ IDX1D4(i+1,j-1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j-1,k-1,s,nx,ny,nz) ] ));
    
    x[MNM] = MNM_BDRY(CMPLX( xr[ IDX1D4(i-1,j  ,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j  ,k-1,s,nx,ny,nz) ] ));
    x[NNM] = NNM_BDRY(CMPLX( xr[ IDX1D4(i  ,j  ,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j  ,k-1,s,nx,ny,nz) ] ));
    x[PNM] = PNM_BDRY(CMPLX( xr[ IDX1D4(i+1,j  ,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j  ,k-1,s,nx,ny,nz) ] ));
    
    x[MPM] = MPM_BDRY(CMPLX( xr[ IDX1D4(i-1,j+1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j+1,k-1,s,nx,ny,nz) ] ));
    x[NPM] = NPM_BDRY(CMPLX( xr[ IDX1D4(i  ,j+1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j+1,k-1,s,nx,ny,nz) ] ));
    x[PPM] = PPM_BDRY(CMPLX( xr[ IDX1D4(i+1,j+1,k-1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j+1,k-1,s,nx,ny,nz) ] ));
    
    x[MMN] = MMN_BDRY(CMPLX( xr[ IDX1D4(i-1,j-1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j-1,k  ,s,nx,ny,nz) ] ));
    x[NMN] = NMN_BDRY(CMPLX( xr[ IDX1D4(i  ,j-1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j-1,k  ,s,nx,ny,nz) ] ));
    x[PMN] = PMN_BDRY(CMPLX( xr[ IDX1D4(i+1,j-1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j-1,k  ,s,nx,ny,nz) ] ));
    
    x[MNN] = MNN_BDRY(CMPLX( xr[ IDX1D4(i-1,j  ,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j  ,k  ,s,nx,ny,nz) ] ));
    x[NNN] =          CMPLX( xr[ IDX1D4(i  ,j  ,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j  ,k  ,s,nx,ny,nz) ] );
    x[PNN] = PNN_BDRY(CMPLX( xr[ IDX1D4(i+1,j  ,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j  ,k  ,s,nx,ny,nz) ] ));
    
    x[MPN] = MPN_BDRY(CMPLX( xr[ IDX1D4(i-1,j+1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j+1,k  ,s,nx,ny,nz) ] ));
    x[NPN] = NPN_BDRY(CMPLX( xr[ IDX1D4(i  ,j+1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j+1,k  ,s,nx,ny,nz) ] ));
    x[PPN] = PPN_BDRY(CMPLX( xr[ IDX1D4(i+1,j+1,k  ,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j+1,k  ,s,nx,ny,nz) ] ));
    
    x[MMP] = MMP_BDRY(CMPLX( xr[ IDX1D4(i-1,j-1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j-1,k+1,s,nx,ny,nz) ] ));
    x[NMP] = NMP_BDRY(CMPLX( xr[ IDX1D4(i  ,j-1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j-1,k+1,s,nx,ny,nz) ] ));
    x[PMP] = PMP_BDRY(CMPLX( xr[ IDX1D4(i+1,j-1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j-1,k+1,s,nx,ny,nz) ] ));
    
    x[MNP] = MNP_BDRY(CMPLX( xr[ IDX1D4(i-1,j  ,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j  ,k+1,s,nx,ny,nz) ] ));
    x[NNP] = NNP_BDRY(CMPLX( xr[ IDX1D4(i  ,j  ,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j  ,k+1,s,nx,ny,nz) ] ));
    x[PNP] = PNP_BDRY(CMPLX( xr[ IDX1D4(i+1,j  ,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j  ,k+1,s,nx,ny,nz) ] ));
    
    x[MPP] = MPP_BDRY(CMPLX( xr[ IDX1D4(i-1,j+1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i-1,j+1,k+1,s,nx,ny,nz) ] ));
    x[NPP] = NPP_BDRY(CMPLX( xr[ IDX1D4(i  ,j+1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i  ,j+1,k+1,s,nx,ny,nz) ] ));
    x[PPP] = PPP_BDRY(CMPLX( xr[ IDX1D4(i+1,j+1,k+1,s,nx,ny,nz) ], xi[ IDX1D4(i+1,j+1,k+1,s,nx,ny,nz) ] ));
}

inline void nbrhood_update(double complex coef[27],double * yr, double * yi, int i, int j,int k, int nx, int ny, int nz,int s,pml_info p)
{
    MMM_IF( yr[ IDX1D4(i-1,j-1,k-1,s,nx,ny,nz) ] += creal(coef[MMM]); yi[ IDX1D4(i-1,j-1,k-1,s,nx,ny,nz) ] += cimag(coef[MMM]);	)
	
    NMM_IF( yr[ IDX1D4(i  ,j-1,k-1,s,nx,ny,nz) ] += creal(coef[NMM]); yi[ IDX1D4(i  ,j-1,k-1,s,nx,ny,nz) ] += cimag(coef[NMM]); )

    PMM_IF( yr[ IDX1D4(i+1,j-1,k-1,s,nx,ny,nz) ] += creal(coef[PMM]); yi[ IDX1D4(i+1,j-1,k-1,s,nx,ny,nz) ] += cimag(coef[PMM]); )
	
    MNM_IF( yr[ IDX1D4(i-1,j  ,k-1,s,nx,ny,nz) ] += creal(coef[MNM]); yi[ IDX1D4(i-1,j  ,k-1,s,nx,ny,nz) ] += cimag(coef[MNM]); )

    NNM_IF( yr[ IDX1D4(i  ,j  ,k-1,s,nx,ny,nz) ] += creal(coef[NNM]); yi[ IDX1D4(i  ,j  ,k-1,s,nx,ny,nz) ] += cimag(coef[NNM]); )

    PNM_IF( yr[ IDX1D4(i+1,j  ,k-1,s,nx,ny,nz) ] += creal(coef[PNM]); yi[ IDX1D4(i+1,j  ,k-1,s,nx,ny,nz) ] += cimag(coef[PNM]); )
	
    MPM_IF( yr[ IDX1D4(i-1,j+1,k-1,s,nx,ny,nz) ] += creal(coef[MPM]); yi[ IDX1D4(i-1,j+1,k-1,s,nx,ny,nz) ] += cimag(coef[MPM]); )

    NPM_IF( yr[ IDX1D4(i  ,j+1,k-1,s,nx,ny,nz) ] += creal(coef[NPM]); yi[ IDX1D4(i  ,j+1,k-1,s,nx,ny,nz) ] += cimag(coef[NPM]); )
	
   PPM_IF( yr[ IDX1D4(i+1,j+1,k-1,s,nx,ny,nz) ] += creal(coef[PPM]); yi[ IDX1D4(i+1,j+1,k-1,s,nx,ny,nz) ] += cimag(coef[PPM]); )

   MMN_IF( yr[ IDX1D4(i-1,j-1,k  ,s,nx,ny,nz) ] += creal(coef[MMN]); yi[ IDX1D4(i-1,j-1,k  ,s,nx,ny,nz) ] += cimag(coef[MMN]); )

   NMN_IF( yr[ IDX1D4(i  ,j-1,k  ,s,nx,ny,nz) ] += creal(coef[NMN]); yi[ IDX1D4(i  ,j-1,k  ,s,nx,ny,nz) ] += cimag(coef[NMN]); )

   PMN_IF( yr[ IDX1D4(i+1,j-1,k  ,s,nx,ny,nz) ] += creal(coef[PMN]); yi[ IDX1D4(i+1,j-1,k  ,s,nx,ny,nz) ] += cimag(coef[PMN]); )

   MNN_IF( yr[ IDX1D4(i-1,j  ,k  ,s,nx,ny,nz) ] += creal(coef[MNN]); yi[ IDX1D4(i-1,j  ,k  ,s,nx,ny,nz) ] += cimag(coef[MNN]); )

	   yr[ IDX1D4(i  ,j  ,k  ,s,nx,ny,nz) ] += creal(coef[NNN]); yi[ IDX1D4(i  ,j  ,k  ,s,nx,ny,nz) ] += cimag(coef[NNN]);

  PNN_IF( yr[ IDX1D4(i+1,j  ,k  ,s,nx,ny,nz) ] += creal(coef[PNN]); yi[ IDX1D4(i+1,j  ,k  ,s,nx,ny,nz) ] += cimag(coef[PNN]); )

  MPN_IF( yr[ IDX1D4(i-1,j+1,k  ,s,nx,ny,nz) ] += creal(coef[MPN]); yi[ IDX1D4(i-1,j+1,k  ,s,nx,ny,nz) ] += cimag(coef[MPN]); )

  NPN_IF( yr[ IDX1D4(i  ,j+1,k  ,s,nx,ny,nz) ] += creal(coef[NPN]); yi[ IDX1D4(i  ,j+1,k  ,s,nx,ny,nz) ] += cimag(coef[NPN]); )

  PPN_IF( yr[ IDX1D4(i+1,j+1,k  ,s,nx,ny,nz) ] += creal(coef[PPN]); yi[ IDX1D4(i+1,j+1,k  ,s,nx,ny,nz) ] += cimag(coef[PPN]); )

  MMP_IF( yr[ IDX1D4(i-1,j-1,k+1,s,nx,ny,nz) ] += creal(coef[MMP]); yi[ IDX1D4(i-1,j-1,k+1,s,nx,ny,nz) ] += cimag(coef[MMP]); )

  NMP_IF( yr[ IDX1D4(i  ,j-1,k+1,s,nx,ny,nz) ] += creal(coef[NMP]); yi[ IDX1D4(i  ,j-1,k+1,s,nx,ny,nz) ] += cimag(coef[NMP]); )

  PMP_IF( yr[ IDX1D4(i+1,j-1,k+1,s,nx,ny,nz) ] += creal(coef[PMP]); yi[ IDX1D4(i+1,j-1,k+1,s,nx,ny,nz) ] += cimag(coef[PMP]); )

  MNP_IF(
    yr[ IDX1D4(i-1,j  ,k+1,s,nx,ny,nz) ] += creal(coef[MNP]);
    yi[ IDX1D4(i-1,j  ,k+1,s,nx,ny,nz) ] += cimag(coef[MNP]);
		)
	    NNP_IF(
    yr[ IDX1D4(i  ,j  ,k+1,s,nx,ny,nz) ] += creal(coef[NNP]);
    yi[ IDX1D4(i  ,j  ,k+1,s,nx,ny,nz) ] += cimag(coef[NNP]);
		)
	    PNP_IF(
    yr[ IDX1D4(i+1,j  ,k+1,s,nx,ny,nz) ] += creal(coef[PNP]);
    yi[ IDX1D4(i+1,j  ,k+1,s,nx,ny,nz) ] += cimag(coef[PNP]);
		)
	    MPP_IF(
    yr[ IDX1D4(i-1,j+1,k+1,s,nx,ny,nz) ] += creal(coef[MPP]);
    yi[ IDX1D4(i-1,j+1,k+1,s,nx,ny,nz) ] += cimag(coef[MPP]);
		)
	    NPP_IF(
    yr[ IDX1D4(i  ,j+1,k+1,s,nx,ny,nz) ] += creal(coef[NPP]);
    yi[ IDX1D4(i  ,j+1,k+1,s,nx,ny,nz) ] += cimag(coef[NPP]);
		)
	    PPP_IF(
    yr[ IDX1D4(i+1,j+1,k+1,s,nx,ny,nz) ] += creal(coef[PPP]);
    yi[ IDX1D4(i+1,j+1,k+1,s,nx,ny,nz) ] += cimag(coef[PPP]);       
		)
}


inline void load_nbrhoodr( double * x, const double * xr, int i, int j, int k, int nx, int ny, int nz, pml_info p )
{
    x[MMM] = MMM_BDRY(xr[ IDX1D4(i-1,j-1,k-1,0,nx,ny,nz) ]);
    x[NMM] = NMM_BDRY(xr[ IDX1D4(i  ,j-1,k-1,0,nx,ny,nz) ]);
    x[PMM] = PMM_BDRY(xr[ IDX1D4(i+1,j-1,k-1,0,nx,ny,nz) ]);
    x[MNM] = MNM_BDRY(xr[ IDX1D4(i-1,j  ,k-1,0,nx,ny,nz) ]);
    x[NNM] = NNM_BDRY(xr[ IDX1D4(i  ,j  ,k-1,0,nx,ny,nz) ]);
    x[PNM] = PNM_BDRY(xr[ IDX1D4(i+1,j  ,k-1,0,nx,ny,nz) ]);
    x[MPM] = MPM_BDRY(xr[ IDX1D4(i-1,j+1,k-1,0,nx,ny,nz) ]);
    x[NPM] = NPM_BDRY(xr[ IDX1D4(i  ,j+1,k-1,0,nx,ny,nz) ]);
    x[PPM] = PPM_BDRY(xr[ IDX1D4(i+1,j+1,k-1,0,nx,ny,nz) ]);
    x[MMN] = MMN_BDRY(xr[ IDX1D4(i-1,j-1,k  ,0,nx,ny,nz) ]);
    x[NMN] = NMN_BDRY(xr[ IDX1D4(i  ,j-1,k  ,0,nx,ny,nz) ]);
    x[PMN] = PMN_BDRY(xr[ IDX1D4(i+1,j-1,k  ,0,nx,ny,nz) ]);
    x[MNN] = MNN_BDRY(xr[ IDX1D4(i-1,j  ,k  ,0,nx,ny,nz) ]);
    x[NNN] =          xr[ IDX1D4(i  ,j  ,k  ,0,nx,ny,nz) ];
    x[PNN] = PNN_BDRY(xr[ IDX1D4(i+1,j  ,k  ,0,nx,ny,nz) ]);
    x[MPN] = MPN_BDRY(xr[ IDX1D4(i-1,j+1,k  ,0,nx,ny,nz) ]);
    x[NPN] = NPN_BDRY(xr[ IDX1D4(i  ,j+1,k  ,0,nx,ny,nz) ]);
    x[PPN] = PPN_BDRY(xr[ IDX1D4(i+1,j+1,k  ,0,nx,ny,nz) ]);
    x[MMP] = MMP_BDRY(xr[ IDX1D4(i-1,j-1,k+1,0,nx,ny,nz) ]);
    x[NMP] = NMP_BDRY(xr[ IDX1D4(i  ,j-1,k+1,0,nx,ny,nz) ]);
    x[PMP] = PMP_BDRY(xr[ IDX1D4(i+1,j-1,k+1,0,nx,ny,nz) ]);
    x[MNP] = MNP_BDRY(xr[ IDX1D4(i-1,j  ,k+1,0,nx,ny,nz) ]);
    x[NNP] = NNP_BDRY(xr[ IDX1D4(i  ,j  ,k+1,0,nx,ny,nz) ]);
    x[PNP] = PNP_BDRY(xr[ IDX1D4(i+1,j  ,k+1,0,nx,ny,nz) ]);
    x[MPP] = MPP_BDRY(xr[ IDX1D4(i-1,j+1,k+1,0,nx,ny,nz) ]);
    x[NPP] = NPP_BDRY(xr[ IDX1D4(i  ,j+1,k+1,0,nx,ny,nz) ]);
    x[PPP] = PPP_BDRY(xr[ IDX1D4(i+1,j+1,k+1,0,nx,ny,nz) ]);
}

inline void load_wn_nbrhood(wn_type wn_window[27],const double * wnr,const double * wni, int i, int j, int k, int nx, int ny, int nz, pml_info p)
{
#ifdef WNCMPLX
    load_nbrhoodc(wn_window,wnr,wni,i,j,k,nx,ny,nz,0,p);
#else
    load_nbrhoodr(wn_window,wnr,i,j,k,nx,ny,nz,p);
#endif 
}

/******************************** 
   
   Coefficient computation functions

 ********************************/

inline coef_consts compute_coef_consts(const double * h)
{
/*
  Compute constants associated to stencil computations
 */
    double hx = h[0]; double hy = h[1]; double hz = h[2]; 
    double hx2 = hx*hx, hy2 = hy*hy, hz2 = hz*hz;
    double hxy = hx2 +hy2; 
    double hxz = hx2 + hz2;
    double hyz = hy2 + hz2;
    double hxyz = hx2 + hy2 + hz2;
    
    double W3A = (W3)*3/(4*hxyz);
    double W3A_2 = 2*W3A;

    double wn_coef  = -(W1 + 3*W2 + 16*W3A*hxyz/3 + WM1-1);
    double wn_xcoef = (W1/hx2 + W2/hx2 + W2/hxz + W2/hxy + 8*W3A);
    double wn_ycoef = (W1/hy2 + W2/hy2 + W2/hyz + W2/hxy + 8*W3A);
    double wn_zcoef = (W1/hz2 + W2/hz2 + W2/hxz + W2/hyz + 8*W3A);

    double pmlx_coef = -(W1/hx2 + W2/hx2 + W2/hxz + W2/hxy + 8*W3A);    
    double pmly_coef = -(W1/hy2 + W2/hy2 + W2/hyz + W2/hxy + 8*W3A);
    double pmlz_coef = -(W1/hz2 + W2/hz2 + W2/hxz + W2/hyz + 8*W3A);

    double xz_coef   = W2/(2*hxz);
    double xy_coef   = W2/(2*hxy);       
    double yz_coef   = W2/(2*hyz);
    
    coef_consts c = { .W3A_2 = W3A_2, .wn_coef = wn_coef, 
		      .wn_xcoef = wn_xcoef, .wn_ycoef = wn_ycoef, .wn_zcoef = wn_zcoef, 
		      .pmlx_coef = pmlx_coef, .pmly_coef = pmly_coef, .pmlz_coef = pmlz_coef, .xz_coef = xz_coef, .xy_coef = xy_coef, .yz_coef = yz_coef };
    return c;
}




inline void get_coefs(double complex coef[27], const wn_type wn_window[27], const coef_consts c, const pml_info p, const pml_adj_info padj)
{
#ifndef ADJ
    /* Compute coefficients - forward mode */
    coef[MMM] = MMM_BDRY(- WM4*wn_window[MMM] + non_deriv_mode*(-c.W3A_2*( p.pmlx_lo + p.pmly_lo + p.pmlz_lo )));
    coef[NMM] = NMM_BDRY(- WM3*wn_window[NMM] + non_deriv_mode*(-c.yz_coef * (p.pmlz_lo + p.pmly_lo) + c.W3A_2*p.pmlx));
    coef[PMM] = PMM_BDRY(- WM4*wn_window[PMM] + non_deriv_mode*(-c.W3A_2*( p.pmlx_hi + p.pmly_lo + p.pmlz_lo )));
    
    coef[MNM] = MNM_BDRY( - WM3*wn_window[MNM] + non_deriv_mode*(-c.xz_coef * (p.pmlz_lo + p.pmlx_lo) + c.W3A_2*p.pmly));
    coef[NNM] = NNM_BDRY(- WM2*wn_window[NNM] + non_deriv_mode*(c.pmlz_coef*p.pmlz_lo + c.yz_coef*p.pmly + c.xz_coef*p.pmlx));
    coef[PNM] = PNM_BDRY(- WM3*wn_window[PNM] + non_deriv_mode*(-c.xz_coef * (p.pmlz_lo + p.pmlx_hi) + c.W3A_2*p.pmly));
    
    coef[MPM] = MPM_BDRY(- WM4*wn_window[MPM] + non_deriv_mode*(-c.W3A_2*( p.pmlx_lo + p.pmly_hi + p.pmlz_lo )));
    coef[NPM] = NPM_BDRY(- WM3*wn_window[NPM] + non_deriv_mode*(-c.yz_coef * (p.pmly_hi + p.pmlz_lo) + c.W3A_2*p.pmlx));
    coef[PPM] = PPM_BDRY(- WM4*wn_window[PPM] + non_deriv_mode*(-c.W3A_2*( p.pmlx_hi + p.pmly_hi + p.pmlz_lo )));		
    
    coef[MMN] = MMN_BDRY(- WM3*wn_window[MMN] + non_deriv_mode*(-c.xy_coef * (p.pmlx_lo + p.pmly_lo) + c.W3A_2*p.pmlz));
    coef[NMN] = NMN_BDRY(- WM2*wn_window[NMN] + non_deriv_mode*(c.pmly_coef*p.pmly_lo + c.yz_coef*p.pmlz + c.xy_coef*p.pmlx));	
    coef[PMN] = PMN_BDRY(- WM3*wn_window[PMN] + non_deriv_mode*(-c.xy_coef * (p.pmlx_hi + p.pmly_lo) + c.W3A_2*p.pmlz));	
    
    coef[MNN] = MNN_BDRY(- WM2*wn_window[MNN]+non_deriv_mode*(c.pmlx_coef*p.pmlx_lo + c.xz_coef*p.pmlz + c.xy_coef*p.pmly));
    coef[NNN] = c.wn_coef*wn_window[NNN] + non_deriv_mode*(c.wn_xcoef*p.pmlx + c.wn_ycoef*p.pmly + c.wn_zcoef*p.pmlz);
    coef[PNN] = PNN_BDRY(- WM2*wn_window[PNN] + non_deriv_mode*(c.pmlx_coef*p.pmlx_hi + c.xz_coef*p.pmlz + c.xy_coef*p.pmly));
    
    coef[MPN] = MPN_BDRY(- WM3*wn_window[MPN] + non_deriv_mode*(-c.xy_coef * (p.pmlx_lo + p.pmly_hi) + c.W3A_2*p.pmlz));		
    coef[NPN] = NPN_BDRY(- WM2*wn_window[NPN] + non_deriv_mode*(c.pmly_coef*p.pmly_hi + c.yz_coef*p.pmlz + c.xy_coef*p.pmlx));
    coef[PPN] = PPN_BDRY(- WM3*wn_window[PPN] + non_deriv_mode*(-c.xy_coef * (p.pmlx_hi + p.pmly_hi) + c.W3A_2*p.pmlz));
    
    coef[MMP] = MMP_BDRY(- WM4*wn_window[MMP] + non_deriv_mode*(-c.W3A_2*( p.pmlx_lo + p.pmly_lo + p.pmlz_hi )));
    coef[NMP] = NMP_BDRY(- WM3*wn_window[NMP] + non_deriv_mode*(-c.yz_coef * (p.pmly_lo + p.pmlz_hi) + c.W3A_2*p.pmlx));
    coef[PMP] = PMP_BDRY(- WM4*wn_window[PMP] + non_deriv_mode*(-c.W3A_2*( p.pmlx_hi + p.pmly_lo + p.pmlz_hi )));
    
    coef[MNP] = MNP_BDRY(- WM3*wn_window[MNP] + non_deriv_mode*(-c.xz_coef * (p.pmlz_hi + p.pmlx_lo) + c.W3A_2*p.pmly));
    coef[NNP] = NNP_BDRY(- WM2*wn_window[NNP] + non_deriv_mode*(c.pmlz_coef*p.pmlz_hi + c.yz_coef*p.pmly + c.xz_coef*p.pmlx));
    coef[PNP] = PNP_BDRY(- WM3*wn_window[PNP] + non_deriv_mode*(-c.xz_coef * (p.pmlz_hi + p.pmlx_hi) + c.W3A_2*p.pmly));
    
    coef[MPP] = MPP_BDRY(- WM4*wn_window[MPP] + non_deriv_mode*(-c.W3A_2*( p.pmlx_lo + p.pmly_hi + p.pmlz_hi )));
    coef[NPP] = NPP_BDRY(- WM3*wn_window[NPP] + non_deriv_mode*(-c.yz_coef * (p.pmlz_hi + p.pmly_hi) + c.W3A_2*p.pmlx));
    coef[PPP] = PPP_BDRY(- WM4*wn_window[PPP] + non_deriv_mode*(-c.W3A_2*( p.pmlx_hi + p.pmly_hi + p.pmlz_hi )));

#else
    /* Compute coefficients - adjoint mode */
    coef[MMM] = - WM4*wn_window[NNN];
    coef[NMM] = - WM3*wn_window[NNN];
    coef[PMM] = - WM4*wn_window[NNN];
    
    coef[MNM] = - WM3*wn_window[NNN];
    coef[NNM] = - WM2*wn_window[NNN];
    coef[PNM] = - WM3*wn_window[NNN];
		
    coef[MPM] = - WM4*wn_window[NNN];
    coef[NPM] = - WM3*wn_window[NNN];
    coef[PPM] = - WM4*wn_window[NNN];		

    coef[MMN] = - WM3*wn_window[NNN];
    coef[NMN] = - WM2*wn_window[NNN];	
    coef[PMN] = - WM3*wn_window[NNN];	

    coef[MNN] = - WM2*wn_window[NNN];
    coef[NNN] = c.wn_coef*wn_window[NNN];
    coef[PNN] = - WM2*wn_window[NNN];
		
    coef[MPN] = - WM3*wn_window[NNN];		
    coef[NPN] = - WM2*wn_window[NNN];
    coef[PPN] = - WM3*wn_window[NNN];
		
    coef[MMP] = - WM4*wn_window[NNN];
    coef[NMP] = - WM3*wn_window[NNN];
    coef[PMP] = - WM4*wn_window[NNN];

    coef[MNP] = - WM3*wn_window[NNN];
    coef[NNP] = - WM2*wn_window[NNN];
    coef[PNP] = - WM3*wn_window[NNN];

    coef[MPP] = - WM4*wn_window[NNN];
    coef[NPP] = - WM3*wn_window[NNN];
    coef[PPP] = - WM4*wn_window[NNN];

#ifndef DERIV
    // (+1,0,0) coef for (i-1,j,k)
    coef[MNN] += c.pmlx_coef*padj.pmlx_hi_window[0] + c.xz_coef*p.pmlz + c.xy_coef*p.pmly;
    // (-1,0,0) coef for (i+1,j,k)
    coef[PNN] += c.pmlx_coef*padj.pmlx_lo_window[2] + c.xz_coef*p.pmlz + c.xy_coef*p.pmly;
    // (0,+1,0) coef for (i,j-1,k)
    coef[NMN] += c.pmly_coef*padj.pmly_hi_window[0] + c.yz_coef*p.pmlz + c.xy_coef*p.pmlx;
    // (0,-1,0) coef for (i,j+1,k)
    coef[NPN] += c.pmly_coef*padj.pmly_lo_window[2] + c.yz_coef*p.pmlz + c.xy_coef*p.pmlx;	
    // (0,0,+1) coef for (i,j,k-1)
    coef[NNM] += c.pmlz_coef*padj.pmlz_hi_window[0] + c.yz_coef*p.pmly + c.xz_coef*p.pmlx;		
    // (0,0,-1) coef for (i,j,k+1)
    coef[NNP] += c.pmlz_coef*padj.pmlz_lo_window[2] + c.yz_coef*p.pmly + c.xz_coef*p.pmlx;		
    // (0,-1,-1) coef for (i,j+1,k+1)
    coef[NPP] += -c.yz_coef * (padj.pmlz_lo_window[2] + padj.pmly_lo_window[2]) + c.W3A_2*p.pmlx;		
    // (0,+1,+1) coef for (i,j-1,k-1)
    coef[NMM] += -c.yz_coef * (padj.pmlz_hi_window[0] + padj.pmly_hi_window[0]) + c.W3A_2*p.pmlx;
    // (0,-1,+1) coef for (i,j+1,j-1)
    coef[NPM] += -c.yz_coef * (padj.pmly_lo_window[2] + padj.pmlz_hi_window[0]) + c.W3A_2*p.pmlx;
    // (0,+1,-1) coef for (i,j-1,j+1)
    coef[NMP] += -c.yz_coef * (padj.pmly_hi_window[0] + padj.pmlz_lo_window[2]) + c.W3A_2*p.pmlx;
    // (-1,0,-1) coef for (i+1,j,k+1)
    coef[PNP] += -c.xz_coef * (padj.pmlz_lo_window[2] + padj.pmlx_lo_window[2]) + c.W3A_2*p.pmly;	
    // (+1,0,+1) coef for (i-1,j,k-1)
    coef[MNM] += -c.xz_coef * (padj.pmlz_hi_window[0] + padj.pmlx_hi_window[0]) + c.W3A_2*p.pmly;				
    // (-1,0,+1) coef for (i+1,j,k-1)
    coef[PNM] += -c.xz_coef * (padj.pmlz_hi_window[0] + padj.pmlx_lo_window[2]) + c.W3A_2*p.pmly;
    // (+1,0,-1) coef for (i-1,j,k+1)
    coef[MNP] += -c.xz_coef * (padj.pmlz_lo_window[2] + padj.pmlx_hi_window[0]) + c.W3A_2*p.pmly;
    // (-1,+1,0) coef for (i+1,j-1,k)
    coef[PMN] += -c.xy_coef * (padj.pmlx_lo_window[2] + padj.pmly_hi_window[0]) + c.W3A_2*p.pmlz;	
    // (+1,-1,0) coef for (i-1,j+1,k)
    coef[MPN] += -c.xy_coef * (padj.pmlx_hi_window[0] + padj.pmly_lo_window[2]) + c.W3A_2*p.pmlz;		
    // (-1,-1,0) coef for (i+1,j+1,k)
    coef[PPN] += -c.xy_coef * (padj.pmlx_lo_window[2] + padj.pmly_lo_window[2]) + c.W3A_2*p.pmlz;		
    // (+1,+1,0) coef for (i-1,j-1,k)
    coef[MMN] += -c.xy_coef * (padj.pmlx_hi_window[0] + padj.pmly_hi_window[0]) + c.W3A_2*p.pmlz;
    // (-1,-1,-1) coef for (i+1,j+1,k+1)
    coef[PPP] += -c.W3A_2*( padj.pmlx_lo_window[2] + padj.pmly_lo_window[2] + padj.pmlz_lo_window[2] );
    // (+1,+1,+1) coef for (i-1,j-1,k-1)
    coef[MMM] += -c.W3A_2*( padj.pmlx_hi_window[0] + padj.pmly_hi_window[0] + padj.pmlz_hi_window[0] );
		
    // (-1,-1,+1) coef for (i+1,j+1,k-1)
    coef[PPM] += -c.W3A_2*( padj.pmlx_lo_window[2] + padj.pmly_lo_window[2] + padj.pmlz_hi_window[0] );		
    // (+1,+1,-1) coef for (i-1,j-1,k+1)
    coef[MMP] += -c.W3A_2*( padj.pmlx_hi_window[0] + padj.pmly_hi_window[0] + padj.pmlz_lo_window[2] );		
    // (-1,+1,-1) coef for (i+1,j-1,k+1)
    coef[PMP] += -c.W3A_2*( padj.pmlx_lo_window[2] + padj.pmly_hi_window[0] + padj.pmlz_lo_window[2] );		
    // (+1,-1,+1) coef for (i-1,j+1,k-1)
    coef[MPM] += -c.W3A_2*( padj.pmlx_hi_window[0] + padj.pmly_lo_window[2] + padj.pmlz_hi_window[0] );		
    // (-1,+1,+1) coef for (i+1,j-1,k-1)
    coef[PMM] += -c.W3A_2*( padj.pmlx_lo_window[2] + padj.pmly_hi_window[0] + padj.pmlz_hi_window[0] );		
    // (+1,-1,-1) coef for (i-1,j+1,k+1)
    coef[MPP] += -c.W3A_2*( padj.pmlx_hi_window[0] + padj.pmly_lo_window[2] + padj.pmlz_lo_window[2] );				    
    // (0,0,0) coef for (i,j,k)
    coef[NNN] += c.wn_xcoef*p.pmlx + c.wn_ycoef*p.pmly + c.wn_zcoef*p.pmlz;

#endif			
#endif
    
    // Boundary handling
    

#ifdef ADJ
    for(int t=0; t<27; t++)
    {
       coef[t] = conj(coef[t]);
    }
#endif


}

void do_Hmvp( const double * wnr, const double * wni, const double * h, const double * n, const double * npml, double *yr, double *yi, const double *xr, const double *xi, int nthreads) {      	 
    int i,j,k,kout,t;
    
    int q;
           
    int nx = (int)n[0]; int ny = (int)n[1]; int nz = (int)n[2];
    int npmlx_lo = (int)npml[0]; int npmlx_hi = (int)npml[1]; 
    int npmly_lo = (int)npml[2]; int npmly_hi = (int)npml[3]; 
    int npmlz_lo = (int)npml[4]; int npmlz_hi = (int)npml[5];

    coef_consts c = compute_coef_consts(h);    

    int pmlz_alloc = 0; int pmly_alloc = 0; int pmlx_alloc = 0;

    double complex coef[27]; 
    double complex x[27]; 
    
    wn_type wn_window[27];

    pml_info p;    
    pml_adj_info padj;
	
    double complex y_out;    
    p.x_hasL = 0; p.x_hasR = 1;
    p.y_hasL = 0; p.y_hasR = 1;
    p.z_hasR = 0; p.z_hasR = 1;
    
#pragma omp parallel for schedule(static) private(coef,x,wn_window,y_out,kout,i,j,k,t,q) firstprivate(p,padj,pmlx_alloc,pmly_alloc,pmlz_alloc) num_threads(nthreads)
    xyzloop_up(
	/* Cache a window of the wavenumber around the current point */	
	load_wn_nbrhood(wn_window,wnr,wni,i,j,k,nx,ny,nz,p);

	/* Get coefficients */
	get_coefs(coef,wn_window, c, p, padj);
	
	/* Cache a window of the wavefield around the current point */
	load_nbrhoodc(x,xr,xi,i,j,k,nx,ny,nz,0,p);

	kout = IDX1D3(i,j,k,nx,ny,nz);
	y_out = 0.0 + 0.0*I;

	for(t=0; t<27; t++){
	    y_out += coef[t] * x[t];
	}
	
	yr[kout] = creal(y_out);
	yi[kout] = cimag(y_out);	    
	)
}
