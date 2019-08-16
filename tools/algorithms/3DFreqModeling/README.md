# 3D Frequency-Domain Modeling Kernel
##  DESCRIPTION
 This package contains the Matlab functions to manipulate the discrete 3D 
 acoustic wave equation operator in frequency domain as well as various quantities 
 depending on solutions to the acoustic wave equation (least squares objective/gradient,
 demigration/migration operators, etc.). This package uses components from tools/algorithms/CommonFreqModeling, as well as 

helm3d_operto27_mvp.m
: Provides an interface to the 27pt stencil Helmholtz matrix vector product (with adjoint + derivative opertions), calls MEX functions

helm3d_operto27_kaczswp.m
: Provides an interface to the 27pt stencil Helmholtz kaczmarz sweep, calls MEX functions

helm3d_7pt_mvp.m
: Provides an interface to the 27pt stencil Helmholtz matrix vector product (with adjoint + derivative opertions)

-- Wrappers to PDEfunc_dist.m / PDEfunc.m
F3d.m
: This function generates 3D acoustic frequency domain data, as well as an associated Jacobian operator. 

opDF3d.m | oppDF3d.m
: SPOT operator for the demigration/migration operator (Jacobian)

opH3D.m | oppH3d.m
: SPOT operator for the full Hessian of FWI

-- Depreciated, explicit matrix code
helmholtz_3d_7p.m - standard 7p matrix stencil

helmholtz_3d.m - 27pt compact stencil of Operto '07

Hmvp.m | Hmvp_MT_mex - C-based code for explicit compressed band-storage matrix vector multiplication



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

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 <http://dx.doi.org/10.1190/1.2759835 [1]> Operto, S., J. Virieux, P. Amestoy, 
 J.-Y. L’Excellent, L. Giraud, and H. B. H. Ali, 2007, 3D finite-difference 
 frequency-domain modeling of visco-acoustic wave propagation using a 
 massively parallel direct solver: A feasibility study: Geophysics, 72, no. 5, 
 SM195–SM211
 
 <http://ethesis.inp-toulouse.fr/archive/00001221/01/pinel.pdf [2]> Pinel, X.,
 2010, A perturbed two-level preconditioner for the solution of 
 three- dimensional heterogeneous Helmholtz problems with applications to 
 geophysics, INPT PhD Thesis
