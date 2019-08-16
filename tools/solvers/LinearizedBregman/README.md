Linearized Bregman Method
##  DESCRIPTION
 This package contains the linearized Bregman method, an iterative solver for L1/L2 constrained problems. 
 The linearized Bregman method solves a version of the basis pursuit (or basis pursuit denoise) problem 
 with a strongly convex objective function
	
	minimize   lambda*||x||_1 + 0.5*||x||_2^2
	subject to A*x = b (or ||A*x-b||_2 < sigma)

For certain lambdas, that are higher than 'some' threshold value, solving the former problem also
solves the basis pursuit (bp denoise) problem.
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
##  DOCUMENTATION
 Documenation to can be accessed by typing 'help lbm' in the matlab terminal.
##  RUNNING
 For solving the basis pursuit (denoise) problem with a strongly convex objective function call

	[x, residual, info] = lbm(A,b,lambda,sigma,x0,options)

The necessary input parameters, as well as the output variables are explained in the function documentation.
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Running applications/demos
 'sparseRecovery.m' shows how to use the function lbm.m to recover a sparse vector from noisy measurements.
 'traceInterpol.m' demonstrates how to use the linearzed Bregman method for seismic trace interpolation.
 using curvelets. Both scripts can be run from the matlab terminal.
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
