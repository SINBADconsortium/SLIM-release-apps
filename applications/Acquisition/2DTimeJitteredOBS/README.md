# 2D ocean-bottom marine acquisition via jittered sampling


## DESCRIPTION

This package is an application of a 2D time-jittered (or blended)
marine acquisition scheme, and a deblending algorithm based on sparse
inversion (via L1 minimization).


## ON-LINE DOCUMENTATION

https://slim.gatech.edu/SoftwareDemos/applications/Acquisition/2DTimeJitteredOBS


## COPYRIGHT

You may use this code only under the conditions and terms of the
license contained in the file LICENSE or COPYING.txt provided with
this source code. If you do not agree to these terms you may not use
this software.


## PREREQUISITES

All prerequisites, except for MATLAB, are provided with the software
release and should be installed before using any of SINBAD's software.


## INSTALLATION

Follow the instructions in the INSTALLATION file (located in the root
directory of this software release) to install all 3-rd party software
(except for MATLAB) and SINBAD's software.


## DOCUMENTATION

Documentation and examples in .html format are included in the ./doc directory.


## RUNNING

Start matlab from this directory or run startup.m to set the correct
paths. The scripts that produced the examples can be found in
./examples/TimeJitAcq_1boat and ./examples/TimeJitAcq_2boats
subdirectories.

1. To run the examples, first download the data from
ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Common/GulfOfSuez by typing
`scons' in the ./data directory

2. Specific instructions for each example can be found in the file README.txt in the corresponding subdirectories.

3. Pre-run results can be downloaded from
ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Acquisition/results by
typing `scons' in the ./results directory

4. To run the algorithm on your own data, see the file
README_templates.txt and the template.m files in the ./scripts
directory.


## SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


## REFERENCES
    
[1] Haneet Wason, and Felix J. Herrmann. Time-jittered ocean bottom
seismic acquisition. Submitted to the SEG (2013).

[2] Haneet Wason, and Felix J. Herrmann. Ocean bottom seismic
acquisition via jittered sampling. EAGE Expanded Abstracts (2013).

[3] Hassan Mansour, Haneet Wason, Tim T.Y. Lin, and Felix
J. Herrmann. Randomized marine acquisition with compressive sampling
matrices. Geophysical Prospecting (2012), vol. 60, 648-662.
    

## NOTE: 

As of October 2012, a correction has been made in the setup
of the time-jittered marine acquisition scheme: speed of the source
vessel is set to the speed during a realistic acquisition scenario,
and is kept constant while the source vessel fires at jittered
instances in time. Hence, the entire setup is in fact a pragmatic
scenario. This was not the case in reference [3], however, the
principles of compressed sensing and the proposed ideas in [3] still
apply after making the correction (see the references cited in [1] and
[2]).

