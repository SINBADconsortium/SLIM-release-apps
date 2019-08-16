# Wavefield Reconstruction Imaging
##  DESCRIPTION
 This application creates images of the interfaces in the medium, given an estimate of the velocity model (background velocity). It has the same goal as other migration methods, but the method is new. It is based on the WRI objective function and is described in [1,2] and an example is shown in [3].
##  ON-LINE DOCUMENTATION
 https://slim.gatech.edu/SoftwareDemos/applications/Imaging/WRimaging/
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
##  DOCUMENTATION
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Downloading data
Download the input model from our ftp server by running the 'Sconstruct' script in in SLIM-release-apps/applications/Imaging/WRimaging/data.
###  Running applications/demos
After downloading the model (or using you own) you can just run the example code (WRImaging_example.m) located in SLIM-release-apps/applications/Imaging/WRimaging/
####  Hardware requirements
* This script was tested using Matlab 2013a with the parallel computing
toolbox.

* Parallelism is achieved by factorizing overdetermined systems (one for each
frequency) in parallel. Each factorization requires about 8 GB. The code will work with any number of factorizations in parallel.
 
* Runtime is about 45 minutes when factorizing 25 overdetermined systems in
parallel on 5 nodes (5 per node). Tested using 2.6GHz Intel processors. 

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES

 [1]<http://dx.doi.org/10.1093/gji/ggt258> 
 Tristan van Leeuwen, Felix J. Herrmann, Geophysical Journal International,2013. Mitigating local minima in full-waveform  inversion by expanding the search space.

 [2]<https://slim.gatech.edu/content/penalty-method-pde-constrained-optimization>
 Tristan van Leeuwen, Felix J. Herrmann. 2013. A penalty method for PDE-constrained optimization.

 [3]<https://slim.gatech.edu/content/examples-penalty-method>
 Bas Peters, Felix J. Herrmann, Tristan van Leeuwen. EAGE, 2014. Wave-equation
based inversion with the penalty method: adjoint-state versus wavefield-reconstruction inversion.
