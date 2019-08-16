# MATLAB interface to running iWAVE
##  DESCRIPTION
This is the toolbox for the MATLAB interface to run iWAVE++. It supports the applications in 'applications/SoftwareDemos/iWAVE'.


##  COPYRIGHT
 You may use this code only under the conditions and terms of the
 license contained in the files LICENSE or COPYING provided with this
 source code. If you do not agree to these terms you may not use this
 software.

##  PREREQUISITES
Using this software requires installing iWAVE software (see INSTALLATION section below). All other prerequisites, except for MATLAB, are provided in the software
release repositories and should be installed as necessary before using
 any of SINBAD's software.

##  INSTALLATION

###  Software in SLIM-release-comp repository
 This application requires installation of extra
 software from SLIM-release-comp repository.

###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary
 components.

###  Installing iWAVE
In order to use this application, you need to install iWAVE distributed with Madagascar version 9520 on your system. (You might want to check with your IT staff if the appropriate version of iWAVE is already installed on your system.) Here is the link that you can use to obtain it from SVN repository using the following command:
 
		svn co -r 9520 https://github.com/ahay/src/trunk/trip iWAVEv9520 
 
Installation of iWAVE depends on the type of MPI environment and compilers. Please, refer to instructions included in obtained source code.

###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.

#### Additional setup for iWAVE
The following elements must be added to shell's PATH environment:

		path_to_iWAVEv9520/rvl/seq/main
		path_to_iWAVEv9520/iwave/acd/main
		path_to_iWAVEv9520/iwave/asg/main
		path_to_iWAVEv9520/iwave/base/main
		path_to_iWAVEv9520/iwave/helm/main
		path_to_iWAVEv9520/iwave/trace/main
		path_to_iWAVEv9520/iwave++/acd++/main
		path_to_iWAVEv9520/iwave++/asg++/main

where path_to_iWAVEv9520 is an absolute path to home of your iWAVE installation.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
