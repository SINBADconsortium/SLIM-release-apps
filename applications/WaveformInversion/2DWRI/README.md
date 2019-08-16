# 2D Wavefield Reconstruction Inversion
##  DESCRIPTION
 This application solves the waveform inversion problem using the Wavefield Reconstruction Inversion method.
##  ON-LINE DOCUMENTATION
https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI/
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
 Follow the instructions in the INSTALLATION file (located in the home directory of this software repository) to install    necessary components.

##  RUNNING

###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Downloading data
Download the input model from our ftp server by running the 'Sconstruct' script in in SLIM-release-apps/applications/WaveformInversion/2DWRI/data.
###  Running applications/demos
After downloading the model you can run the example code (exp_wri_BGcompass_full_offset.m) located in SLIM-release-apps/applications/WaveformInversion/2DWRI/. 'plot_WRI_example_results.m' can be used for plotting results.
####  Hardware requirements
* This script was tested using Matlab 2015b with the parallel computing
toolbox.

* Parallelism is achieved by factorizing overdetermined systems (one for each
frequency) in parallel. Each factorization requires about 15 GB. The example is computed most efficiently with 4 factorizations in parallel, because each frequency batch in the example has 4 frequencies. Additional parallelism can be achieved by using the multi-threaded version of the factorization algorithm, which can be activated in the example script. The total number of cores required is the number of frequencies factored in parallel times the number of threads used per factorization. Factoring each overdetermined system on a different node/core is intrinsically a parallel operation, whereas multi-threaded factorization of a single system is not.
 
* Runtime is about 45 minutes hours when factorizing 4 overdetermined systems in
parallel. Tested using 2.6GHz Intel processors. 

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
