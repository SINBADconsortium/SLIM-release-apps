# 2D Wavefield Reconstruction Inversion with Source Estimation
##  DESCRIPTION
This application implements the source estimation for wavefield reconstruction inversion.
##  ON-LINE DOCUMENTATION
 You can find the on-line documentation of this application by clicking the following link:

  https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI-SourceEstimation/
  
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
 This application does not require installation of extra
 software from SLIM-release-comp repository.
###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary components.
##  RUNNING
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Running applications/demos
 You can run the example code (`Example_source_estimation.m`) located in `SLIM-release-apps/applications/WaveformInversion/2DWRI-SourceEstimation`. After you open an **Matlab**, you should first start the parallel environment by the following command:

 `pool=parpool(3);`

 This command will open a matlab pool with 3 workers. Then you can run the code by the command:

 `Example_source_estimation`.


####  Hardware requirements
* This script was tested using Matlab 2015b with the parallel computing
toolbox.

* Parallelism is achieved by factorizing overdetermined systems (one for each frequency) in parallel. Each factorization requires about 1 GB. The example is computed with 3 factorizations in parallel.

* Runtime is about 1 hour when factorizing 3 overdetermined systems in parallel. Tested using 2.6GHz Intel processors.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 1. https://slim.gatech.edu/content/source-estimation-wavefield-reconstruction-inversion-0 Zhilong Fang, Rongrong Wang, and Felix J. Herrmann, Submitted to Geophysics, 2017, Source estimation for wavefield-reconstruction inversion.
