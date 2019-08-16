#Fast and robust full-waveform inversion without source encoding.
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
This applications is available only in the software release for members of SINBAD consortium.

##DESCRIPTION
This package contains matlab scripts that reproduce examples from [1,2]
    using robust penalties, source-estimation and fast optimization without source-encoding. 
    It also includes a small demo that reproduces the famous Camembert example [3].
    An overview of the results can be found in the ./doc directory and online: <br />
    <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DRobustBatching/>
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
##DOCUMENTATION
Documentation and examples in .html format are included in the ./doc subdirectory.
    Help can also be accessed from Matlab via the `help' command.
    
##RUNNING
Start matlab from this directory or run startup.m to set the correct paths. 
The scripts that produced the examples shown in the documentation can be found in ./scripts.

1. To run the scripts, first download the data and models from
    	<ftp://ftp.slim.gatech.edu/data/users/tristan/2DRobustBatching/data>
       by typing `scons' in the ./data directory
2. Specific instructions for each script can be found in the scripts. 
3. Pre-run results can be downloaded from
	<ftp://ftp.slim.gatech.edu/data/users/tristan/2DRobustBatching/results>
       by typing `scons' in the ./results directory
4. To run the algorithms on your own data, see template.m in the ./scripts directory. 

##NOTES
A part of the Compass Model, developed by the BG group is distributed with this release. 
    Please refer to BG_DISCLAIMER.txt in the ./data directory before downloading the data.
##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##REFERENCES

	[1] T. van Leeuwen and F.J. Herrmann - Fast Waveform inversion without source-encoding, Geophysical Prospecting, submitted
	[2] A.Y. Aravkin, T. van Leeuwen and F.J. Herrmann - Source estimation for frequency-domain FWI with robust penalties, EAGE Expanded abstracts 2012, submitted.
	[3] O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)
