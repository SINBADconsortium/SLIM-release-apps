# Time domain LSRTM with sparsity promotion

##  DESCRIPTION
 This Demo performs least squares RTM with sparsity promotion and source subsampling in the time domain. The optimization problem that we solve with our sparsity promoting LSRTM algorithm is:
	
	minimize   lambda*||C*x||_1 + 0.5*||C*x||_2^2
	subject to J*x = d (or ||J*x-d||_2 < sigma)

C is the 2-dimensional curvelet transform, i.e. the sparsity is enforced in the curvelet domain. J is the linearized modeling operator (Born modeling), d is the observed, linearized data and x is the unknown model perturbation/image. The problem is solved via the linearized Bregman method using a subsampled set of sources in each iteration, which significantly reduces the overall number of PDE solves and makes this algorithm computationally feasible.

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
 components.

##  RUNNING

###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.

###  Download input data
 An example data set for the Marmousi model can be downloaded from the ftp server by running the 'Sconstruct' script.

###  Run data example
 First, open MATLAB and load the startup file from the current folder (where README.md is). Then change directory (in MATLAB) to examples folder.
 The demo file in the examples folder performs LSRTM on the Marmousi model. For an arbitrary model and data set, the user needs to supply an input data set (ideally linearized data) and a smooth background velocity model. Several data and model parameters need to be supplied inside the demo script, as well as the parameters for the LSRTM (e.g. number of iterations and shots per iteration).
 Runtime for the Marmousi example is approximately 12 hours (40 iterations with 8 sources per iteration in parallel).

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
