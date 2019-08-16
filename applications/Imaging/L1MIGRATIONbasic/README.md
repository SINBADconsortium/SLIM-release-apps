# 2D Efficient least-squares imaging with sparsity promotion and compressive sensing.

##DESCRIPTION
The computational cost for solving least-squares imaging problem is prohibitive. 
	It mainly depends on the number of partial-differential equation solves required 
	while inversion. In this work, we borrow ideas from compressive sensing and 
	stochastic optimization to turn the original over-determined least-squares 
	imaging problem into a under-determined problem. With this approach, it can 
	allow us to compute least-squared imaging with a subsampled data, yielding a 
	significant speed-up.

This package is an application for the Efficient least-squares imaging. In 
	which we include a toy model example as well as examples with BG compass model 
	(synthetic model created with constraint from real well-log information).
    An overview of the results can be found in the ./Doc directory.
## ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/Imaging/L1MIGRATIONbasic>
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
	
Before using this package you should make sure your environment is
setup correctly. First run 'matlab\_path.m' which in the home directory
of this package in Matlab to append all the included software 
modules to the search path.
## DOCUMENTATION
Documentation and examples in .html format are included in the ./Doc subdirectory.<br />
Help can also be accessed from Matlab via the `help' command.
## RUNNING
This package contains the examples in [1] which can be found in './Examples/'. <br />

1. To run the examples, first go to './Data/' folder, and run<br />
		'scons' to fetch the data from SLIM ftp
2. Then go to './Examples/' folder, run 'LSM_toy_model_L1_WR.m' using 3 workers in
                parallel pool to test the algorithm;
                or run 'LSM_BGcompass_model2d_L1_WR.m' using 10 workers in parallel pool
                for the result based on BG model.
3. To view the results, please go to './Results/' folder. Then
		you can find saved results here. Since this package using an
		iterative method, you can view the results iteratively, and also
		you can run 'scons' under './Results/' folder to get pre-run results.
		
## NOTES
A part of the Compass Model, developed by the BG group is distributed
with this release.<br />
Please refer to BG_DISCLAIMER.txt in the ./Data directory before
downloading the data.
## SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

## REFERENCES
	[1] Felix J. Herrmann and Xiang Li, “Efficient least-squares imaging with sparsity promotion and compressive sensing”, Geophysical Prospecting, vol. 60, p. 696-712, 2012.
	[2] Felix J. Herrmann and Xiang Li, “Randomized dimensionality reduction for full-waveform inversion”, in EAGE Technical Program Expanded Abstracts, 2010.
