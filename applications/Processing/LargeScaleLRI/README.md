# Large scale, parallel low-rank interpolation using SPGLR
##  DESCRIPTION
 This package illustrates the main features of large scale seismic data completion using SPGLR, the low-rank matrix completion solver based on LR-factorizations of the underlying matrix and the SPGL1 framework. The scripts in these directories are the seismic examples associated to the SPGLR_PAR package, which can be found in /tools/solvers/SPGLR_PAR/.
 
 ##  ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/Processing/LargeScaleLRI/>

##  COPYRIGHT
 You may use this code only under the conditions and terms of the license contained in the files LICENSE or COPYING provided with this source code. If you do not agree to these terms you may not use this software.

##  PREREQUISITES
 All prerequisites, except for MATLAB, are provided in the software  release repositories and should be installed as necessary before using any of SINBAD's software. This particular software package requires the parallel computing toolbox for Matlab, as well as a system configuration that supports parallel workers.

##  INSTALLATION
 Follow the instructions in the INSTALLATION file (located in the home directory of this software repository) to install necessary components.

##  DOCUMENTATION
 Documentation and examples in .html format are included in the ./doc subdirectory. Help can also be accessed from Matlab via the 'help' command.

##  RUNNING
Start Matlab from this directory or run startup.m to set the correct paths. The scripts that produce the examples can be found in ./examples/. To download the example data, run 'scons' from ./data directory. The parameters for each experiment can be found in the ./spgLR_experiments/ directory. You can define your own experiment by copying + renaming one of the experiment files with the name spgLR_experiment##.m, where ## is the unique ID number. 

Example usage:
experiment = 1;
spgLR_bgdata; %will run matrix completion interpolation with paramters specified in spgLR_experiment1.m


###  Preparing shell environment
You must setup your shell environment according to the steps listed in the README located in home directory of the software release.
###  Downloading data
 Run the SConstruct scons scripts in the ./data/ and ./results/ to download the data files as well as a sample results file for a single provided experiment.

###  Running applications/demos
The parameters for each matrix completion experiment are encapsulated in the spgLR_experiment###.m files in the ./spgLR_experiments/ subdirectory. You can create a new experiment with new parameters by cloning one of these files, giving it the name spgLR_experiment###.m with ### replaced by a unique number, and running the spgLR_bgdata.m file in the ./examples subdirectory with the new ### you gave your experiment file as a parameter.

####  Hardware requirements
   This code is meant to be run in parallel using matlab's parallel software toolbox. A multi-core or multi-node system is required to run this code.   
   
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
1. C. Da Silva, F. Herrmann "Large-scale seismic data interpolation in a parallel computing environment." Available at <https://slim.gatech.edu/content/large-scale-seismic-data-interpolation-parallel-computing-environment>.
2. A. Y. Aravkin, et al. "Fast methods for denoising matrix completion formulations, with applications to robust seismic data interpolation", SIAM Journal on Scientific Computing, 2014. Available at <https://slim.gatech.edu/content/fast-methods-denoising-matrix-completion-formulations-applications-robust-seismic-data>.
