##  DESCRIPTION
 This package illustrates the main features of the 3D modeling code
 that is found under /tools/algorithms/3DFreqModeling and tools/solvers/Krylov
 
##  ON-LINE DOCUMENTATION
A very comprehensive formatted online README along with a short execution walkthrough
is available here:<br />
 <https://slim.gatech.edu/SoftwareDemos/applications/Modeling/3DAcousticFreqModeling/>
It is a good idea to refer to the above URL first and only refer to this
'README.txt' and 'README_examples.txt' files for detailed technical information.

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
 Follow the instructions in the INSTALLATION file (located in the
 root directory of this software release) to install all 3-rd party
 software (except for MATLAB) and SINBAD's software.
 
##  DOCUMENTATION
 Documentation and examples in .html format are included in the ./doc subdirectory.  
 Help can also be accessed from Matlab via the `help' command.
 
##  RUNNING
 Start matlab from this directory or run startup.m to set the correct paths. 
 The scripts that produced the examples shown in the documentation can be found 
 in ./scripts and in ./examples.

 To properly run some of the demos and tests, it is necessary to download data
 files from the SLIM FTP server. 
 
 *Installation with SCONS*
 Enter the data directory and use the scons executable:
     
     cd ./data
     scons
     
 *Installation without SCONS*
 Download all the data files found in 
 
     ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Modeling/3DAcousticFreqModeling/data
     
 and place them into ./data folder.
 
 
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
