# Time-domain Gauss-Newton full-waveform inversion for Chevron 2014 benchmark dataset
##  DESCRIPTION
This package contains matlab scripts that generates some inversion result for 
Chevron 2014 SEG workshop benchmark dataset. This is an application of SLIM
2D/3D time-stepping forward modeling kernel[1] which can be found in the tools 
folder of this software release. 
##  ON-LINE DOCUMENTATION
An overview of the results can be found in the ./doc directory and online: <br />
 <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DBasicTimeStepping/>
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
###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary
 components.
##  DOCUMENTATION
Documentation and examples in .html format are included in the ./doc subdirectory.
Help can also be accessed from Matlab via the `help' command.
##  RUNNING
Start matlab from this directory or run startup.m to set the correct paths. 
The scripts that produced the examples shown in the documentation can be found in ./scripts.
###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.
###  Downloading data
please download Chevron 2014 SEG workshop benchmark dataset from 
<https://s3.amazonaws.com/open.source.geoscience/open_data/seg_workshop_fwi_2014/seg_workshop_fwi_2014.html>
You can also find DATA LICENSE AGREEMENT information with the above link. Please
use this dataset according to thire requirement presented in the above link.
###  Running applications/demos
Specific instructions for each script can be found in the scripts. 
####  Hardware requirements
To run the script, you need at least 40 cores. each core needs a least 2GB memory. 
##  NOTES
Please use Chevron 2014 SEG workshop benchmark dataset according to their DATA LICENSE AGREEMENT, we are not
response for any consequence of using that data, if you do NOT follow Chevron license agreement.
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
##  REFERENCES
	[1]. Xiang Li - 2D/3D time-stepping for imaging and inversion, SINBAD Fall consortium talks. 2014

