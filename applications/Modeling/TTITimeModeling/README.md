# 2D/3D acoustic anisotropic time-domain modeling and linearized modeling
##  DESCRIPTION
This package contains basic demonstration for 2D TTI modelling and 3D acoustic modelling with domain decomposition
##  ON-LINE DOCUMENTATION
https://slim.gatech.edu/SoftwareDemos/applications/Modeling/TTITimeModeling/
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
###  Software in SLIM-release-comp repository
 This application requires installation of extra
 software from SLIM-release-comp repository.

	SPOT - object oriented framework for matrix-free linear algebra.
	pSPOT - parallel extension of SPOT.

operators/misc

	opLInterp1D - 1D cubic Lagrange interpolation
	opExtension - Pads input with zeros or constant values
	opSmooth - 1D smoothing by convolution with triangular kernel
	opSaveSnapshot - save history of iterations operator
functions/misc

    grid2odn, odn2grid - convert grid vectors to [origin, increment, size] triplet and vice versa
###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary
 comp[onents.
##  DOCUMENTATION
You can view the documentation for the example at 
https://slim.gatech.edu/SoftwareDemos/applications/Modeling/TTITimeModeling/
##  RUNNING
In order to run the example you only need to go in the scripts directory and directly run one of the two scripts.
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Downloading data
This example doesn't need data.
####  Hardware requirements
In order to run this example you need the parallel toolbox (for domain decomposition). You do not need more than 8Gb of memory as this example is a small 2D and 3D layer model. Both the example runs within 20 minutes.
####  Data adaptation
You can change the basic layer model by any velocity or anisotropy model. If you do so you will need to change 

	model.n,model.d and model.o according to your model
	The acquisition geometry (model.xsrc,....)
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
##  REFERENCES
 1. https://slim.gatech.edu/Publications/Private/TechReport/2015/witte2015TRoam/witte2015TRoam.html 
