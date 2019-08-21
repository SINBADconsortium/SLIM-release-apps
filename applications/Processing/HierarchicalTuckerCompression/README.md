# Large-scale seismic data compression with on-the-fly shots/receivers generation from compressed Hierarchical Tucker parameter


## DESCRIPTION 
This application illustrates the fast shots/recievers extraction directly from Hierarchiacal Tucker parameters once we either compress or interpolate the data in HT format, which can lead to substantial reduction in memory costs in subsequent shot-based processings.

## ON-LINE DOCUMENTATION 
You can find the on-line documentation of this application in the following link:
 
<https://slim.gatech.edu/SoftwareDemos/applications/Processing/HierarchicalTuckerCompression/>


## COPYRIGHT 
You may use this code only under the conditions and terms of the license contained in the files LICENSE or COPYING provided with this source code. If you do not agree to these terms you may not use this software.


## PREREQUISITES 

This application does not require installation of extra packages from SLIM-release-comp repository.


## INSTALLATION 
Follow the instructions in the INSTALLATION file (located in the home directory of this software repository) to install necessary components.

## DOCUMENTATION 
Documentation and examples in .html format are included in the ./doc subdirectory.
Help can also be accessed from Matlab via the `help' command.

## RUNNING 
From the home of this application directory, you can follow several steps to get the results:
 
1. Change directory to ./data and Fetch data from the FTP server by calling 'scons' inside this directory. Return to application's home when done.

1. Start Matlab in home directory of this application to run startup.m and to set the correct MATLAB paths. 

3. Change directory to ./examples and run the 'BG_3D.m' script to obtain the results (and you do not need parallel matlab workers in this example).

4. The results would have been already stored in the ./results after finishing running the 'BG_3D.m' script, and this script would plot all those results on your screen (self-contained)

## SUPPORT 
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


## REFERENCES 
1. Y. Zhang, C. Da Silva, R. Kumar and F. J. Herrmann, "Massive 3D seismic data compression and inversion with hierarchical Tucker", SEG Technical Program Expanded Abstracts, 2017. Available at <https://slim.gatech.edu/content/massive-3d-seismic-data-compression-and-inversion-hierarchical-tucker>
2. C. Da Silva, Y. Zhang, R. Kumar and F. J. Herrmann, "Applications of low-rank compressed seismic data to full waveform inversion and extended image volumes", Submitted to Geophysics, 2018. Available at <https://slim.gatech.edu/content/applications-low-rank-compressed-seismic-data-full-waveform-inversion-and-extended-image>
