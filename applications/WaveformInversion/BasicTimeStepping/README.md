# Time-domain 2D FWI with TTI anisotropy
##  DESCRIPTION
This demo is a basic example of time-domain full-waveform inversion on the 2D BG compass model. We use gradient descent with line-search and we show the result for acoustic isotropic inversion and acoustic anisotropic inversion.
##  ON-LINE DOCUMENTATION
https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/BasicTimeStepping/
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
 components.
##  DOCUMENTATION
You can view the documentation for the example at 
https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/BasicTimeStepping/FWI_TTI_2D.html
##  RUNNING
To use this software
	add the data folder to the path
	run the script FWI_TTI_2D.m
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Downloading data
	run the Sonstruct script in the data directory to get the necessary dataset
####  Hardware requirements
 This example is serial but contains an spmd block. If you use Matlab 2015b a parpool will be open and the source experiments will be done in parallel but the script will run properly. If you don't have a  parallel toolbox, you need to deactivate the automatic parpool in the preferences of Matlab.
 This is a simple 2D example and it should run in 20 minutes on a single CPU and you need at least 8Gb of RAM to store the wavefields.
####  Data adaptation
If you wish to run the example on a different model :
	replace the velocity and anisotropy files by yours with the appropriate reading commands
	modify model.n, model.d and model.o according to your model.
	Set the source and receiver (model.xsrc,...) position to your choice.
	The true data will be generated for the true velocity model. If you don not have a true velocity model you need to replace the lines 
	
		[mm,modelm,~,anim]=Setup_CFL(m,model,[],ani);
		q=sp_RickerWavelet(modelm.f0,1/modelm.f0,modelm.dt,modelm.T);
		dataT=Gen_data(mm,modelm,q,[],anim);

by your own data.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 1. https://slim.gatech.edu/Publications/Private/TechReport/2015/witte2015TRoam/witte2015TRoam.html 
