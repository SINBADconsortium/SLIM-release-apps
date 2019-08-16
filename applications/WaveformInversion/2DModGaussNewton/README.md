#2D modified Gauss-Newton full-waveform inversion.
##DESCRIPTION
This package is an application for modified gauss-newton (GN) 			
	full-waveform inversion, which based on the ideas from 
	compressive-sensing and stochastic optimization, where the GN updates 
	are computed from random subsets of the data via curvelet-domain 
	sparsity-promotion. There are two different subset sampling strategies 
	are considered in this package: randomized source encoding, and drawing 
	sequential shots firing at random source locations from marine data 
	with missing near and far offsets. In both cases, we obtain excellent 
	inversion results compared to conventional methods at reduced 
	computational costs. There is also a small example based on Camembert
	example [3] which can allow users to test the algorithm in a short time 
    An overview of the results can be found in the ./doc directory.
##ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DModGaussNewton>
##COPYRIGHT
You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.
##PREREQUISITES
All prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
##INSTALLATION
Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.
	
Before using this package you should make sure your environment is
	setup correctly. First run 'matlab_path.m' which in the home directory
	of this package in Matlab to append all the included software 
	modules to the search path.
##DOCUMENTATION
Documentation and examples in .html format are included in the ./doc subdirectory.
  	Help can also be accessed from Matlab via the `help' command.
##RUNNING
This package contains the examples in [1,2] which can be found in
'./examples/'.

1. To run the examples, first go to './data/' folder, and run 
		'scons' to fetch the data from SLIM ftp
2. Then go to './examples/example\_camembert/' folder, run 'MGNFWI\_camenbert.m' to test
		the algorithm; or go to './examples/example_BG/' folder, run 'MGNFWI_BG.m' for the
		result based on BG model.
3. To view the results, please go to './results/' folder. Then
		you can find saved results here. Since this package using an
		iterative method, you can view the results iteratively, and also
		you can run 'scons' under './results/' folder to get pre-run results.
		
##NOTES
A part of the Compass Model, developed by the BG group is distributed 	
	with this release. 
 	Please refer to BG_DISCLAIMER.txt in the ./data directory before 
	downloading the data.
##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##REFERENCES

	[1] Felix J. Herrmann, Xiang Li, Aleksandr Y. Aravkin, and Tristan van Leeuwen, “A modified, sparsity promoting, Gauss-Newton algorithm for seismic waveform inversion”, in Proc. SPIE, 2011, vol. 2011.
	[2] Xiang Li, Aleksandr Y. Aravkin, Tristan van Leeuwen, and Felix J. Herrmann, “Fast randomized full-waveform inversion with compressive sensing”. 2011. Geophysics, accepted.
	[3] O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)
