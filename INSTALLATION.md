# INSTALLATION instructions
## 1 INTRODUCTION
 We provided simple installation scripts for each software package for
 your convenience; however, we cannot guarantee that they will be
 sufficient to successfully complete installation. Please, contact
 SLIM's developers (see SUPPORT section) with any questions or problems
 encountered during installation.
 
 The 3-rd party software packages included with SLIM software release
 might require additional prerequisites depending on the version of the
 operating system and type and/or completeness of its installation.
 This installation was tested and completes successfully on full
 installation of x86_64 Linux with operating systems CentOS 5.5 or
 later (equivalent to RedHat 5.5 or later).
 
 If any step of the installation fails for the 3-rd party software, it
 might be more effective to consult first your IT support if some of
 the prerequisites are missing or the included web site for specific
 instructions suitable for your operating system.
## 2 INSTALLATION STRATEGIES
 Here are the potential strategies for installing our software in
 single-/multi-user environment.
 
 Note, that each user must have a personal (or at the very least has to
 have write permissions) copy of SLIM-release-apps in order to to run
 any applications, since each application is configured to look for
 data and create directories inside of its own directory.
### 2.1 Complete user-owned installation
 Personal installation for single user, both SLIM-release-comp and
 SLIM-release-apps.
 
 User installs both SLIM-release-comp and SLIM-release-apps in personal
 directories.
### 2.2 Multi-user SLIM-release-comp installation
 Multi-user installation of SLIM-release-comp and personal installation
 of SLIM-release-apps. 
 
 Designated person instals SLIM-release-comp in common directory and
 user installs personal SLIM-release-apps.
### 2.3 Full multi-user installation
 Multi-user installation of SLIM-release-comp, multi-user installation
 of SLIM-release-apps, and personal copy of SLIM-release-apps
 
 Designated person instals both SLIM-release-comp and SLIM-release-apps
 in common directories, and user obtains personal copy of
 SLIM-release-apps (no installation needed).
## 3 DOWNLOADING
 In terminal, change directory to the location where you want to
 install the software and execute the following git command:
 
 	git clone git@github.com:SINBADconsortium/SLIM-release-apps.git
 
 and the cloned software will be in SLIM-release-apps sub-directory.
## 4 GIT BRANCHES
 SLIM is using master branch to develop/add software to repository. To
 avoid using the software that is not complete or fully tested, you
 might want to use use `stable` branch. You may switch to `stable`
 branch using the following:
 
 	git checkout stable
## 5 SHELL ENVIRONMENT
 You must configure your shell environment before you can proceed with
 installation.
### 5.1 Prepare environment.* script
#### 5.1.1 Make a copy of appropriate environment.* script template
 Open terminal window, change directory to the home of
 SLIM-release-apps and do either of the following:
##### 5.1.1.1 in bash-like shell execute
 	cp environment.sh.template environment.sh
##### 5.1.1.2 in csh-like shell execute
 	cp environment.csh.template environment.csh
#### 5.1.2 Edit and modify your copy of the  environment.* script
 Edit the file created above in an editor and make the following
 changes  accordingly to chosen installation strategy.
 
 Note: if you skipped the installation of SLIM-release-comp, leave
 SLIM_COMP environment as it is in any of the options below.
##### 5.1.2.1 Complete user-owned installation
 In personal SLIM-release-apps/environment.{sh,csh}:
 
 - SLIM_COMP points to absolute path of user's personal
 SLIM-release-comp  
 - SLIM_APPS points to absolute path of user's personal
 SLIM-release-apps  
 - SLIM_APPS_RUNS points to $SLIM_APPS  
 
##### 5.1.2.2 Multi-user SLIM-release-comp installation
 In personal SLIM-release-apps/environment.{sh,csh}:
 
 - SLIM_COMP points to absolute path of common SLIM-release-comp  
 - SLIM_APPS points to absolute path of user's personal
 SLIM-release-apps  
 - SLIM_APPS_RUNS points to $SLIM_APPS  
 
##### 5.1.2.3 Full multi-user installation
 In common SLIM-release-apps/environment.{sh,csh}:
 
 - SLIM_COMP points to absolute path of common SLIM-release-comp  
 - SLIM_APPS points to absolute path of common SLIM-release-apps  
 - SLIM_APPS_RUNS points to $SLIM_APPS  
 
 In personal SLIM-release-apps/environment.{sh,csh}:
 
 - SLIM_COMP points to absolute path of common SLIM-release-comp  
 - SLIM_APPS points to absolute path of common SLIM-release-apps  
 - SLIM_APPS_RUNS points to absolute path of user's personal
 SLIM-release-apps  
 
### 5.2 Importing shell environment
 In terminal window, in the home of SLIM-release-apps do either of the
 following:
#### 5.2.1 in bash-like shell execute
 	. environment.sh
#### 5.2.2 in csh-like shell execute
 	source environment.csh
#### 5.3 Loading the environment automatically (optional and for users only)
 Add the following to your default-shell's startup script to make the
 permanent change to the environment. You will not need then to source
 environment.sh/csh manually.
##### 5.3.1 for bash-like shell add
 	. path_to-SLIM-release-apps/environment.sh
##### 5.3.2 for csh-like shell add
 	source path_to-SLIM-release-apps/environment.csh
### 5.4 Testing the environment
 Once configured, you can check if the SLIM_COMP , SLIM_APPS, and
 SLIM_APPS_RUNS  environments are set correctly and verify if MATLAB
 executables are in the PATH using the following:
#### 5.4.1 in bash-like shell execute
 	test_env4slim.sh
#### 5.4.2 in csh-like shell execute
 	test_env4slim.csh
## 6 MATLAB startup for batch jobs
 Users who intend to run their jobs in non-interactive batch mode must
 add those extra steps. Note, it is not necessary to do those
 modification while only installing the software for other users.
### 6.1 Load environment automatically
 Make sure to follow the steps in the section above to automatically
 source environment.* from shell's startup scripts.
### 6.2 Add startup_slim.m script to those ~/matlab directory
  First. create ~/matlab directory if it does not exists.
 
 Execute the following command in terminal:
 
 	cp $SLIM_APPS/skel/startup_SLIM.m ~/matlab
### 6.3 Call startup_SLIM script at the end of ~/matlab/startup.m
  First, create ~/matlab/startup.m if it does not exist.
 
 Add the following line at the end of ~/matlab/startup.m:
 
 	startup_SLIM
## 7 INSTALLATION
 The installation relies on correctly set shell environment (see SHELL
 ENVIRONMENT section) and should be executed in the terminal window
 only after that particular terminal session have the environment
 configured. Before starting/continuing with the installation, make
 sure that your environment was setup according to the instruction
 above.
### 7.1 MATLAB version R2014a or later (http://www.mathworks.com/)
 We cannot provide here the installation for MATLAB software due to
 MATLAB license limitation. However, please, make sure to have version
 R2014a installed at the very least. The test_env4slim.* (mentioned
 above in SHELL environment section) will check the if necessary MATLAB
 commands are present, and warn you if version of MATLAB is too old.
### 7.2 Compile all MATLAB's mex binaries
 Running installation script does not require much time, and typically
 takes a couple of seconds. The script prints to the screen the
 information about current stage and the name of the output log. This
 information will be needed if you encounter any errors during
 compilation.
 
 To compile MATLAB's mex binaries in SLIM's MATLAB toolboxes execute
 the following in the terminal:
 
 	install_MEX
 
 If you encounter problems during installation (FATAL ERRORS), please,
 remember to attach the log files from the failing stage of
 installation while contacting us for support.
## 8 TESTING
 To complete the installation, one should run an extra script that will
 test MATLAB mex binaries. Contrary to installation scripts, the
 testing script will output its progress to the terminal. Please, check
 the output for potential errors.The comments in the subsections below
 explain the meaning of each test.
 
 To run test suite execute the following from the terminal:
 
 	install_test4matlab
 
 Note, that the testing suite will stop for a few seconds between each
 test.
 
 Contact SLIM's developers if any of those tests fail.
### 8.1 Serial tests
#### 8.1.1 SPOT
 This test runs a series of SPOT unit tests and will report results to
 the terminal. All tests should report `pass` in the 3-rd column.
#### 8.1.2 SPGL1
 We did not devise yet a proper unit test for SPGL1, so SPGL1 part of
 the test executes suite of demos. Some of the demos are expected to
 fail. However the test should neither produce any MATLAB errors nor
 stop during execution.
#### 8.1.3 Quasi-Newton
 This test runs a single Quasi-Newton solver test. It should report
 `PASSED`.
#### 8.1.4 HTOpt
 The test runs a few serial tests in Hierarchical Tucker Optimization
 toolbox. It should report `PASSED`.
#### 8.1.5 LinearizedBregman
 The test for Linearized Bregman method. It should report `PASSED`.
#### 8.1.6 2D Curvelet transform.
 This test runs a dot-test for 2D curvelet interface, and only if
 SLIM-release-comp is installed. It should report `PASSED`.
#### 8.1.7 3D Curvelet transform.
 This test runs a dot-test for 3D curvelet interface, and only if
 SLIM-release-comp is installed. It should report `PASSED`.
#### 8.1.8 MATLAB for RSF
 This a short test to check MATLAB's RSF interface, and only if
 SLIM-release-comp is installed.  It should report `PASSED`.
#### 8.1.9 Tests in tools/algorithms/3DFreqModeling
 Series of tests for algorithm in tools/algorithms/3DFreqModeling. It
 should report `PASSED`.
#### 8.1.10 Tests in tools/algorithms/TimeModeling
 Series of tests for algorithm in tools/algorithms/2DFreqModeling. It
 should report `PASSED`.
#### 8.1.11 Tests in tools/operators/misc
 Unit tests for miscellaneous operators in directory
 tools/operators/misc
### 8.2 Parallel tests
 Tests in this category require parallel MATLAB license. All tests use
 default MATLAB's `parallel pool` configuration with 2 workers.
#### 8.2.1 pSPOT
 This test runs a series of pSPOT unit tests and will report results to
 the terminal. It will utilize 2 workers to check parallel
 functionality. All tests should report `pass` in the 3-rd column.
 
 Note, this test will fail if you do not have working (not installed of
 lacking the license) at the least MATLAB's Parallel Toolbox, or your
 worker configuration is constrained for minimum number of workers
 bigger then 2. Then the test will report that it is not able to open
 properly parallel session, so called 'parallel pool'.
 
#### 8.2.2 HTOpt
 The test runs a few parallel tests in Hierarchical Tucker Optimization
 toolbox. It should report `PASSED`.
#### 8.2.3 Tests in tools/algorithms/2DFreqModeling
 Series of tests for algorithm in tools/algorithms/2DFreqModeling. It
 should report `PASSED`.
#### 8.2.4 Tests in tools/algorithms/3DFreqModeling including tools/solvers/Krylov
 Series of tests for functions in tools/algorithms/3DFreqModeling/
 including tools/solvers/Krylov.  It should report `PASSED`.
#### 8.2.5 Tests in tools/algorithms/TimeModeling
 Series of tests for algorithm in tools/algorithms/2DFreqModeling. It
 should report `PASSED`.
## 9 SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
