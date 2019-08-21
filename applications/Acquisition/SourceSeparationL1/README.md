# Source separation for towed-streamer marine data via sparsity promotion


##  DESCRIPTION

This package contains a MATLAB implementation of a 2-D over/under
blended marine acquisition scheme, and a deblending (or source
separation) algorithm based on sparsity-promoting inversion in the
curvelet domain using L1 minimization.


##  ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/Acquisition/SourceSeparationL1>


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

This application requires installation of extra software from
SLIM-release-comp repository.

###  Software in SLIM-release-apps (this) repository

Follow the instructions in the INSTALLATION file (located in the home
directory of this software repository) to install necessary
components.


##  DOCUMENTATION

Documentation and examples in .html format are included in the ./doc
directory.


##  RUNNING

Start matlab from this directory or run startup.m to set the correct
paths. The scripts that produced the examples can be found in
./examples.

###  Preparing shell environment

You must setup your shell environment according to the steps listed in
the README located in home directory of the software release.

###  Downloading data

To run the examples, first download the data from
ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Acquisition/SourceSeparationL1/Original,
by typing `scons' in the ./data directory

###  Running applications/demos

Specific instructions for the example can be found in the README.md
file in the corresponding subdirectory.

Pre-run results can be downloaded from
ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Acquisition/SourceSeparationL1/Results,
by typing `scons' in the ./results directory

#### Data adaptation
    
To run the algorithm on your own data, see the instructions in the
README.md file in the ./examples directory.


##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


##  REFERENCES

[1] Rajiv Kumar, Haneet Wason, and Felix J. Herrmann [2015]. "Source
separation for simultaneous towed-streamer marine acquisition --- a
compressed sensing approach", Geophysics, vol. 80, pp. WD73-WD88.
Available at
<https://slim.gatech.edu/content/source-separation-simultaneous-towed-streamer-marine-acquisition-–-compressed-sensing>.

