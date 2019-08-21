# 3D Frequency-Domain Modeling Kernel
##  DESCRIPTION
 This package contains the Matlab functions to manipulate the discrete 3D 
 acoustic wave equation operator in frequency domain as well as various quantities 
 depending on solutions to the acoustic wave equation (least squares objective/gradient,
 demigration/migration operators, etc.).
    
helmholtz_3d.m
: Computes the discrete Helmholtz operator using a staggered grid second order
  finite differences scheme. The output is an array with  dimension 
  (27,nx*ny*nz), or in other words, it stores the discrete operator in a band-
  storage format rather than as a sparse matrix. Each column of this array 
  represents a point in the physical domain.
  This stencil is discussed in [1] and it aims at minimizing the 
  dispersion error, attaining stability with as little as 4 points per 
  wavelength in homogeneous test-cases. We usually recommend the use of 6 points
  per wavelength for heterogeneous media.

helmholtz_3d_derivative.m
: Computes the tensor ∂H/∂m, where H is the Helmholtz matrix discretized using 
  helmholtz_3d.m routine, and m is the model. This is required, for instance,
  to compute the gradient of the misfit function efficiently in the adjoint-
  state formulation, instead of computing the full Jacobian. Although this is a
  tensor, it can be stored exactly as H itself, in band-storage format. Using
  this object might be a bit trickier than expected. Please, refer to the full
  documentation, the demos or Matlab's Help on how to use this tensor.
  
Hmvp.m
: Matrix vector multiplication with the output of helmholtz_3d. It calls a mex
  function to perform the matrix-vector product itself in an efficient and 
  multi-threaded environment.

create_pml.m
: Function used exclusively to compute the damping coefficients for the PML 
  layer. This function is called internally by helmholtz_3d and never needs to 
  be called manually by the user.

H2sparse.m
: Converts the output of helmholtz_3d to a standard Matlab's sparse matrix 
  format. Specially useful if you intent to use direct solvers.

Htransp.m
: Transposes the matrix computed using helmholtz_3d (or helmholtz_3d_7p), 
  preserving the band-storage format.
  
helmholtz_3d_7p.m (beta)
: Very similar to helmholtz_3d.m, except that this generates a discrete 
  operator using the standard second order 7 points stencil with the PML; 
  please refer to [2] Appendix A for more details on this stencil. 
  This generates a  matrix with 7 diagonals only (as opposed to 27) but 
  requires more points per  wavelength. In the literature it is common the use 
  of 10 points per wavelength, but 12 points is a more conservative and also 
  popular choice. 
  
 For the sake of compatibility, the legacy version of the 3D modeling code is 
 still maintained, but will be removed in the future (to be announced). The
 list of functions in the legacy version is:
 
    F3d_old     - modeling operator
    DF3d        - Jacobian of modeling operator
    A_Helm3D    - Construct 3D Helmholtz operator as sparse matrix
    G_Helm3D    - Jacobian of A_Helm3D
    oppDF3d     - pSPOT wrapper for DF
    pF3d        - parallel version of modeling operator
    pR_Helm3D.m - parallel version of constructing 3D Helmholtz operator
    R_Helm3D.m  - Constructing 3D Helmholtz operator only storing the relative index and value
    pRtransp.m  - parallel version for the transpose of the Helmholtz operator
    (p)CARP(B)CG - (parallel) (Block) CARP-CG iterative method to solve the Helmholtz equation
    R2mat, mat2R, Rtransp - utility functions for creating and manipulating band-storage matrices
    sweepR_mex  - Kaczmarz sweeps on band-storage matrices, used by (p)CARP(B)CG.
    test_*        - unit tests
    Hmvp_mex    - matrix-vector products with a band-storage matrix

PDEfunc3D.m 
: This function computes various quantities depending on solutions of the Helmholtz equation in 3D such as forward modeling, migration/demigration, gauss-newton hessian, hessian products, least-squares objectives + gradients. This is a serial code that is run on each worker in a parallel environment. The user can specify options as to what solver to use to solve the Helmholtz equations (default: CARPCR). 

PDEfunc3D_dist.m
: This function distributes the full computation of the output for PDEfunc3D over sources/frequencies using the Parallel Matlab Toolbox. 

opBandStorage.m
: This SPOT operator is used to encapsulate the Helmholtz matrix in diagonal storage format. It implements matrix-vector products by calling the appropriate mex functions. Given the proper parameters, it also knows how to divide itself by calling a specified linear solver (with possible preconditioners). 

####Important Remarks (For Developers):

1. This package allows free surface; this is achieved by setting the PML length 
   in the top of the domain to zero. See the example for more details.
2. the old routines generates a (n_x*n_y*n_z, 27); this code inverts this to the
   more efficient (27, n_x*n_y*n_z).
3. the Legacy version uses slowness squared. The new version allows the user to 
   choose between meters per second or slowness squared (this is emphasized in 
   the comments inside the code itself).
4. The legacy code permutes the z direction to the first dimension. This code 
   keeps z direction as the third dimension.

##  ON-LINE DOCUMENTATION
 The online documentation along with some demos and tests can be found at
 <https://slim.gatech.edu/SoftwareDemos/applications/Modeling/3DAcousticFreqModeling/>

##  COPYRIGHT
 You may use this code only under the conditions and terms of the
 license contained in the files LICENSE or COPYING provided with this
 source code. If you do not agree to these terms you may not use this
 software.

##  PREREQUISITES
 All prerequisites, except for MATLAB, are provided in the software
 release repositories and should be installed as necessary before using
 any of SINBAD's software.

##  INSTALLATION
###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary
 comp[onents.

##  RUNNING
 Please, read the online documentation at
 <https://slim.gatech.edu/SoftwareDemos/applications/Modeling/3DAcousticFreqModeling/>
 for examples on how to use these functions. They are readily available from
 Matlab. You can also use "help command" for more information on the
 interface of every function. 

## NOTES
This package does not support parallel computation of the discrete operator as 
of now. It is possible, however, to compute the discrete operator first and then 
distribute it - see the demos for an example on how to do that. We are still
working on a full parallel version of this package.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 <https://doi.org/10.1190/1.2759835> Operto, S., J. Virieux, P. Amestoy, 
 J.-Y. L’Excellent, L. Giraud, and H. B. H. Ali, 2007, 3D finite-difference 
 frequency-domain modeling of visco-acoustic wave propagation using a 
 massively parallel direct solver: A feasibility study: Geophysics, 72, no. 5, 
 SM195–SM211
 
 <http://oatao.univ-toulouse.fr/7217/> Pinel, X.,
 2010, A perturbed two-level preconditioner for the solution of 
 three- dimensional heterogeneous Helmholtz problems with applications to 
 geophysics, INPT PhD Thesis
