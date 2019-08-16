# Time-jittered blended marine acquisition on non-uniform spatial grids

##  DESCRIPTION

This package contains a MATLAB implementation of a 2-D time-jittered
blended marine acquisition scheme on non-uniform spatial (source)
grid, and a deblending algorithm based on sparse inversion via
L1-minimization incorporating the non-equispaced fast discrete
curvelet transform.


##  ON-LINE DOCUMENTATION

https://slim.gatech.edu/SoftwareDemos/applications/Acquisition/TimeJitteredOBS_OffTheGrid


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

Examples in .html format are included in the ./doc directory.


##  RUNNING

Start matlab from this directory or run startup.m to set the correct
paths. The scripts that produced the examples can be found in
./examples.

###  Preparing shell environment

You must setup your shell environment according to the steps listed in
the README located in home directory of the software release.

###  Downloading data

To run the examples, first download the data from
ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Common/GulfOfSuez by typing
`scons' in the ./data directory

###  Running applications/demos

Specific instructions to run the examples can be found in the README.md
file in the ./examples directory.

####  Data adaptation

To run the deblending algorithm on your own data, see the instructions
in the README.md file in the ./examples directory.


##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


##  REFERENCES

[1] Haneet Wason, and Felix J. Herrmann [2013]. "Time-jittered ocean
bottom seismic acquisition." SEG Technical Program Expanded Abstracts,
vol. 32, pp. 1-6. Available at
https://slim.gatech.edu/content/time-jittered-ocean-bottom-seismic-acquisition.

[2] Gilles Hennenfent, Lloyd Fenelon, and Felix J. Herrmann
[2010]. "Nonequispaced curvelet transform for seismic data
reconstruction: a sparsity-promoting approach." Geophysics, vol. 75,
pp. WB203-WB210. Available at
https://slim.gatech.edu/content/nonequispaced-curvelet-transform-seismic-data-reconstruction-sparsity-promoting-approach.

[3] Haneet Wason, Felix Oghenekohwo, and Felix J. Herrmann
[2015]. "Compressed sensing in 4-D marine – recovery of dense
time-lapse data from subsampled data without repetition." Presented at
the 77th EAGE Conference and Exhibition, doi:
10.3997/2214-4609.201413088. Available at
https://slim.gatech.edu/content/compressed-sensing-4-d-marine–-recovery-dense-time-lapse-data-subsampled-data-without-repeti.

