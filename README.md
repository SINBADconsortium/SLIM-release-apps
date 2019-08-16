# Public legacy release of MATLAB-software repository for SLIM software release to SINBAD consortium
## 1 DESCRIPTION
 This is the repository for SLIM's software release for SINBAD
 consortium members. The release contains applications and tools
 (including 3-rd party prerequisites) necessary to demonstrate and use
 algorithms developed by SLIM's researchers. SLIM's software release is
 organized in two repositories:
### 1.1 Repository SLIM-release-apps - this repository
 Repository at [SLIM-release-apps]
 (https://github.com/SINBADconsortium/SLIM-release-apps) contains core
 of SLIM's software, i.e. all applications, algorithms, tools, and
 utilities. The software in this repository requires minimal
 installation, but some of the applications depend on installation of
 SLIM-release-comp listed below.
 
 Note that each user must have a private copy of SLIM-release-apps in
 order to to run any applications, since an application is configured
 to look for data and create directories inside of application's
 directory. See INSTALLATION file for possible installation strategies.
### 1.2 Repository SLIM-release-comp
 Repository at [SLIM-release-comp]
 (https://github.com/SINBADconsortium/SLIM-release-comp) containing
 extra 3rd-part software for multi-user installation - not needed for
 some applications from repository SLIM-release-apps. The installation
 of software from this repository may be shared by multiple users and
 may require lengthy installation.  See INSTALLATION file for possible
 installation strategies.
## 2 COPYRIGHT
 You may use this code only under the conditions and terms of the
 license contained in the file LICENSE provided with this source code.
 If you do not agree to these terms you may not use this software.
## 3 PREREQUISITES
 All prerequisites, except for MATLAB, are provided with the software
 release and should be installed as necessary before using SLIM's
 software. Installation and using SLIM's software requires MATLAB
 version R2014a or later from [MathWorks] (http://www.mathworks.com/)
 including Parallel Computing Toolbox (check file MATLAB-details in
 this repository for more information). In addition, the installation
 requires GNU gcc/g++ compilers that come with the operating system.
 Included in this repository are:
### 3.1 Software packages in SLIM-release-comp repository
 Those packages are listed in README of SLIM-release-comp repository.
### 3.2 Software packages in SLIM-release-apps (this) repository
#### 3.2.1 3rd-party software packages
 [NFFT] (https://www-user.tu-chemnitz.de/~potts/nfft/) - native MATLAB
 implementation  
#### 3.2.2 SLIM's software packages
 SLIM's version of [SPGL1] (http://www.cs.ubc.ca/labs/scl/spgl1/)  
 SLIM's version of [SPOT] (http://www.cs.ubc.ca/labs/scl/spot/)  
 [pSPOT] (https://github.com/slimgroup/pSPOT/)  
## 4 INSTALLATION
 Follow the instructions in the INSTALLATION file in this repository to
 install core of SLIM's software and create necessary scripts to
 configure your environment. If installation of 3rd-part software from
 SLIM-release-comp is required, follow the instructions in INSTALLATION
 file of SLIM-release-comp before you install software in this
 repository. You will need to have MATLAB installed and added to your
 PATH environment before you can proceed. If you encounter any problems
 during the installation, please, let us know. See SUPPORT section at
 the end for contact information.
## 5 SHELL ENVIRONMENT
 You must configure your shell environment before you can run any
 applications. The environment.* scripts, created during installation,
 will set SLIM_COMP (and if needed and available), SLIM_APPS, and
 SLIM_APPS_RUNS environments and add necessary executable and library
 paths to your shell. You will not be able to run applications if
 SLIM_COMPS /SLIM_APPS/SLIM_APPS_RUNS  variables are not set correctly
 or MATLAB executables are missing from your shell's PATH environment.
### 5.1 Importing shell environment
 In order to import shell environment you will have to source one of
 the scripts in the home of SLIM-release-apps. 
#### 5.1.1 From the home of SLIM-release-apps
 In the terminal window and once per terminal session, change directory
 to the home of SLIM-release-apps and do either of the following:
##### 5.1.1.1 in bash-like shell execute
 	. environment.sh
##### 5.1.1.2 in csh-like shell execute
 	source environment.csh
#### 5.1.2 From another location
 If you want to source either of those scripts from the location other
 then the home of SLIM-release-apps (like from your shell's startup
 scripts). Then, in the terminal window (and only once per terminal
 session) do either of the following from any location:
##### 5.1.2.1 in bash-like shell execute
 	. path_to-SLIM-release-apps/environment.sh
##### 5.1.2.2 in csh-like shell execute
 	source path_to-SLIM-release-apps/environment.csh
##### 5.1.2.3 optional for interactive jobs and obligatory for batch submission
 Add the appropriate one of the above to your default-shell's startup
 script to make the permanent change to the environment. You will not
 need then to source environment.sh/csh manually.
### 5.2 Testing your environment
 Once configured, you can check if the SLIM_COMP , SLIM_APPS, and
 SLIM_APPS_RUNS environments are set correctly and verify if MATLAB
 executables are in the PATH using:
#### 5.2.1 in bash-like shell execute
 	test_env4slim.sh
#### 5.2.2 in csh-like shell execute
 	test_env4slim.csh
## 6 MATLAB startup for batch jobs
 Users who intend to run their jobs in non-interactive batch mode must
 add those extra steps:
### 6.1 Add startup_slim.m script to those ~/matlab directory
  First. create ~/matlab directory if it does not exists.
 
 Execute the following command in terminal:
 
 	cp $SLIM_APPS/skel/startup_SLIM.m ~/matlab
### 6.2 Call startup_SLIM script at the end of ~/matlab/startup.m
  First, create ~/matlab/startup.m if it does not exist.
 
 Add the following line at the end of ~/matlab/startup.m:
 
 	startup_SLIM
## 7 RUNNING APPLICATIONS
 Please, see the specific instructions that are included in the README
 file in each application's directories. The application directories
 are organized by topics and are located inside of applications/
 directory in the root of unpacked SLIM's software release. See SLIM
 software's components sections in this README for the current list of
 applications.
 
 Note, that you might need to download the data before running an
 application. The relevant instructions can be found in application's
 README file.
### 7.1 Typical execution sequence
 Running any application (after the input file(s) are downloaded)
 involves typically the following steps:
 
 - open Matlab  
 - start parallel pool with appropriate number of workers, if needed  
 - change directory to the main directory of desired application  
 - run `startup` if `startup.m` file exists in this directory  
 - change directory to the directory of scrips/examples if such exists 
 - run `startup` if `startup.m` file exists in this directory  
 - run desired script  
 
 The above steps will ensure that the toolboxes necessary for the
 application are loaded and the input data files are found.
 
### 7.2 `RunApplication` helper function
 A provided MATLAB function `RunApplication` helps running our
 applications. It executes the steps above (except for starting the
 parallel pool) and ensures that: 1) appropriate toolboxes are added
 before application is being used, and 2) that the application is
 executed in proper location. For more information how to use this
 function type in MATLAB:
 
 	help RunApplication
 
#### 7.2.1 Using `RunApplication` in interactive mode
 You will need to create parallel pool with appropriate number of
 workers before executing `RunApplication`.
#### 7.2.2 Using `RunApplication` in batch (non-interactive) mode
 Function `RunApplication` also allows to easily submit our
 applications to run in the non-interactive batch mode using the
 following syntax:
 
 	batch(@RunApplication, 0, {x1,..., xn}, *other_batch_options...*)
 
 where  `{x1,..., xn}` are the same arguments as used for
 `RunApplication` in the interactive mode, but enclosed in the cell
 array. For more information about batch command, and its options, type
 in MATLAB:
 
 	help batch
 
## 8 DOCUMENTATION
 For the information about SLIM's software packages, please, check the
 README files included with each package. Especially, check the README
 files and the applications to see how to execute and customize the
 applications to change parameters and/or use them with other input
 data. On-line documentation  for all applications can be found at
 [slim.gatech.edu]
 (https://slim.gatech.edu/software-documentation).
## 9 SLIM software's components
 The SLIM's software is divided into to main components: applications
 and tools. Tools are mostly for general use or are typically problem
 independent; e.g., algorithms and solvers are considered as tools.
 Applications illustrate how those different tools can be used to solve
 specific problems.
### 9.1 Applications
#### 9.1.1 applications/Acquisition/2DTimeJitteredOBS
 2D ocean-bottom marine acquisition via jittered sampling
#### 9.1.2 applications/Acquisition/2DTimeJitteredOBS-LR
 Rank minimization based source separation in time-jittered marine
 acquisiotion
#### 9.1.3 applications/Acquisition/SourceSeparationLowRankHSS
 Source separation via SVD-free rank minimization in the hierarchical
 semi-separable representation
#### 9.1.4 applications/Acquisition/TimeJitteredOBS_OffTheGrid
 Time-jittered blended marine acquisition on non-uniform grids
#### 9.1.5 applications/Acquisition/TimeLapseJRM
  Joint recovery method for time-lapse seismic data
#### 9.1.7 applications/Imaging/L1MIGRATIONwSRM
 Fast imaging with surface-related multiples by sparse inversion
#### 9.1.9 applications/Imaging/TimeDomainLSRTM
 Time domain LSRTM with sparsity promotion
#### 9.1.10 applications/Imaging/WRimaging
 Wavefield reconstruction imaging
#### 9.1.11 applications/Modeling/2DAcousticFreqModeling
 Tutorial for 2D Frequency-domain acoustic modelling and imaging
#### 9.1.13 applications/Modeling/TTITimeModeling
 Time-domain 2D/3D modeling and linearized modeling
#### 9.1.14 applications/deprecated/Modeling/3DAcousticFreqModeling
 3D Frequency-Domain Modeling Kernel
#### 9.1.15 applications/Processing/HierarchicalTuckerOptimization
 Missing receiver interpolation of 3D frequency slices using
 Hierarchical Tucker optimization
#### 9.1.16 applications/Processing/LargeScaleLRI
 Large-scale seismic data interpolation using SVD-free low-rank matrix
 factorization.
#### 9.1.17 applications/Processing/LowRankInterpolationAndDenoising
 Seismic data regularization, interpolation, and denoising using
 factorization based low-rank optimization
#### 9.1.18 applications/Processing/SparsityPromotingDenoising
 Sparsity-promoting denoising of seismic data
#### 9.1.19 applications/Processing/WeightedL1Interpolation
 Seismic trace interpolation using weighted one-norm minimization
#### 9.1.20 applications/WavefieldSeparation/RobustEPSI
 Robust Estimation of Primaries by Sparse Inversion (via L1
 minimization)
#### 9.1.21 applications/WaveformInversion/2DBasic
 2D Basic Acoustic Full Waveform Inversion
#### 9.1.24 applications/WaveformInversion/2DRobustBatching
 Fast and robust 2D full-waveform inversion without source encoding.
#### 9.1.25 applications/WaveformInversion/2DWRI
 2D Wavefield Reconstruction Inversion
#### 9.1.26 applications/WaveformInversion/2DWRI-TVconstrained
 Total Variation Regularized Wavefield Reconstruction Inversion
#### 9.1.27 applications/WaveformInversion/3DBasic
 3D FWI with an Acoustic Helmholtz Modeling Kernel
#### 9.1.28 applications/WaveformInversion/3DBatching
 3D acoustic full-waveform inversion
#### 9.1.29 applications/WaveformInversion/3DParallelBatching
 Parallel 3D frequency domain full waveform inversion
#### 9.1.30 applications/WaveformInversion/BasicTimeStepping
 Time-domain 2D FWI with TTI anisotropy
#### 9.1.31 applications/WaveformInversion/ConstrainedFWI
 Constrained FWI
#### 9.1.32 applications/WaveformInversion/TimeLapseFWI
 Full Waveform Inversion for time-lapse seismic data
#### 9.1.33 applications/SoftwareDemos/2DSMII
 Scripts to reproduce the examples of the paper by T. van Leeuwen, 2012.
 ["A parallel matrix-free framework for frequency-domain seismic
 modeling, imaging and inversion."]
 (https://slim.gatech.edu/node/27141)
 
#### 9.1.34 applications/SoftwareDemos/iWAVE
 Examples for using iWAVE interface for different applications
### 9.2 Tools in SLIM-release-apps
#### 9.2.1 tools/algorithms/CommonFreqModeling
 Common tools for frequency-domain acoustic modelling
#### 9.2.2 tools/algorithms/2DFreqModeling
 2D Frequency-domain acoustic modelling
#### 9.2.3 tools/algorithms/3DFreqModeling
 3D Frequency-domain acoustic modelling
#### 9.2.4 tools/algorithms/AdaptiveSparseRecovery
 Adaptive sparse recovery
#### 9.2.5 tools/algorithms/LowRankMinimization
 The algorithm for seismic data regularization, interpolation, and
 denoising using factorization based low-rank optimization
#### 9.2.6 tools/algorithms/REPSI
 Robust Estimation of Primaries by Sparse Inversion.
#### 9.2.7 tools/algorithms/TimeModeling
 Time-domain 2D/3D acoustic modeling kernel
#### 9.2.8 tools/algorithms/WRI
 Wavefield reconstruction inversion
#### 9.2.9 tools/solvers/GenSPGL1
 Generalized SPGL1
#### 9.2.10 tools/solvers/HTOpt
 Hierarchical Tucker Optimization
#### 9.2.11 tools/solvers/Krylov
 Krylov solvers for band-storage operators
#### 9.2.12 tools/solvers/LinearizedBregman
 Linearized Bregman method for L1 constrained optimization.
#### 9.2.13 tools/solvers/Multigrid
 Multigrid preconditioners
#### 9.2.14 tools/solvers/QuasiNewton
 Quasi Newton optimization method: L-BFGS with weak Wolfe linesearch
#### 9.2.15 tools/solvers/SPGL1-SLIM
 SLIMs' version of SPGL1 1.6 - a solver for large-scale sparse
 reconstruction. Adapted for parallel computations.
#### 9.2.16 tools/solvers/SPGLR_PAR
 Parallel SPGLR - large-scale matrix completion
#### 9.2.17 tools/utilities/SPOT-SLIM
 SLIM's version of SPOT - a linear-operator toolbox for MATLAB. Adapted
 for needs of pSPOT.
#### 9.2.18 tools/utilities/pSPOT
 Parallel extensions to SPOT.
#### 9.2.19 tools/utilities/iWAVE
 MATLAB interface to running iWAVE
### 9.3 Tools in SLIM-release-comps
#### 9.3.1 tools/transforms/CurveLab-2.1.2-SLIM
 SLIM's version of CurveLab 2.1.2 - curvelet transform.
## 10 SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
