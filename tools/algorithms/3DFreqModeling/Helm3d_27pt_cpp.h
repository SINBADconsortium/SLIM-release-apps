#include "math.h"
#include <thread>
#include <vector>

using namespace std;

/**
   Enum definitions
 */
enum Mult_Mode {
    FORW_MULT,
    ADJ_MULT
};

enum Deriv_Mode {
    DERIV_MODE,
    NO_DERIV_MODE
};

enum Wavenum_Cmplx {
    WN_IS_CMPLX,
    WN_IS_REAL
};
    
template <class T,Mult_Mode m, Deriv_Mode d, Wavenum_Cmplx c>
void do_Hmvp( const T * wnr, const T * wni, const T * h, const T * n, const T * npml, T * yr, T *yi, const T * xr, const T * xi, int zmin, int zmax);


template <class T,Mult_Mode m, Deriv_Mode d, Wavenum_Cmplx c>
void do_Hmvp_mt( const T * wnr, const T * wni, const T * h, const T * n, const T * npml, T * yr, T * yi, const T * xr, const T * xi, int n_threads){
    vector<thread> th(n_threads);
    int zmin, zmax;    
    
    int dz = (int)(n[2])/n_threads;
   
    for (int i=0; i<n_threads; i++){
        zmin = i*dz;
        if (i<n_threads-1)
            zmax = (i+1)*dz;
        else
            zmax = n[2];
        th[i] = thread(do_Hmvp<T,m,d,c>,wnr,wni,h,n,npml,yr,yi,xr,xi,zmin,zmax);
    }

    for (auto &t : th) {
        t.join();
    }


    
}
