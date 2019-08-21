# Uncertainty quantification for 2D Wavefield Reconstruction Inversion
##  DESCRIPTION
This application implements the uncertainty quantification for wavefield reconstruction inversion.
##  ON-LINE DOCUMENTATION
 You can find the on-line documentation of this application by clicking the following link:

  <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/UQforWRI/>

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

### Download the data and model
To run the scripts, first download the data and models from<br />
       <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/UQforWRI/data/>
      by typing `scons` in the ./data directory

You can also download the result of the script from <br />
       <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/UQforWRI/result/>
       by typing `scons` in the ./result directory

###  Running applications/demos
 You can run the example code (`Example_UQ.m`) located in `SLIM-release-apps/applications/WaveformInversion/UQforWRI`. After you open an **Matlab**, you should first start the parallel environment by the following command:

 `pool=parpool(3);`

 This command will open a matlab pool with 3 workers. Then you can run the code by the command:

 `Example_UQ`.


####  Hardware requirements
* This script was tested using Matlab 2015b with the parallel computing
toolbox.

* Parallelism is achieved by factorizing overdetermined systems (one for each frequency) in parallel. Each factorization requires about 1 GB. The example is computed with 3 factorizations in parallel.

* Runtime is about 6 hour when factorizing 3 overdetermined systems in parallel. Tested using 2.6GHz Intel processors.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 1. <https://slim.gatech.edu/content/uncertainty-quantification-inverse-problems-weak-pde-constraints> Zhilong Fang, Rongrong Wang, and Felix J. Herrmann, Submitted to Geophysics, 2017, Source estimation for wavefield-reconstruction inversion.
