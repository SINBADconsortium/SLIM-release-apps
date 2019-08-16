# Rank minimization based source-separation in time-jittered 2D ocean-bottom marine acquisition.


## DESCRIPTION
This package is an application of a 2D time-jittered (or blended) marine 
     acquisition scheme, and a deblending algorithm based on rank-minimization.
The acquisition setup for generating the blended data is same as proposed by Haneet Wason.
However, we offer an efficient deblending algorithm using rank-minimization based techniques.


## ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/Acquisition/2DTimeJitteredOBS-LR>
## COPYRIGHT
You may use this code only under the conditions and terms of the
     license contained in the file LICENSE or COPYING.txt provided with
     this source code. If you do not agree to these terms you may not
     use this software.


## PREREQUISITES
All prerequisites, except for MATLAB, are provided with the
     software release and should be installed before using any of
     SINBAD's software.


## INSTALLATION
Follow the instructions in the INSTALLATION file (located in the
     root directory of this software release) to install all 3-rd party
     software (except for MATLAB) and SINBAD's software.


## DOCUMENTATION
Documentation and examples in .html format are included in the ./doc directory.


## RUNNING
Start matlab from this directory or run startup.m file from this directory 
to set the correct paths. The scripts that produced the example can be found in ./examples. 

1. To run the example, first download the data from<br />
    	 <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Common/GulfOfSuez><br />
        by typing `scons' in the ./data directory

2. Specific instructions for the example can be found in the file README.md in the 
        ./examples directory.

3. Pre-run results can be downloaded from<br />
        <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Acquisition/2DTimeJitteredOBS-LR/results><br />
        by typing `scons' in the ./results directory

4. To run the algorithm on your own data, see the file README.md in the 
        ./examples directory.


## SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


## REFERENCES
    [1] Rajiv Kumar, Haneet Wason, and Felix J. Herrmann. Time-jittered marine acquisition: 
        low-rank v/s sparsity. Submitted to the EAGE (2015).

    [2] Haneet Wason, and Felix J. Herrmann. Time-jittered ocean bottom seismic acquisition.
        Submitted to the SEG (2013).
    

