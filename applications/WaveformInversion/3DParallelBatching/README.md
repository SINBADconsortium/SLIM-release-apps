#Parallel 3D frequency domain full waveform inversion
This applications is available only in the software release for members of SINBAD consortium.

#DESCRIPTION
This package contains matlab scripts that demonstrate the 3D Parallell FWI framework implemented in 
    A unified 2D/3D software environment for large scale time-harmonic full waveform inversion, Curt Da Silva and Felix Herrmann, 2016.<br />
    An overview of the results can be found in the ./doc directory and online:<br />
    <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/3DParallelBatching/>
#COPYRIGHT
You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.
#PREREQUISITES
All prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
#INSTALLATION
Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.
#DOCUMENTATION
Documentation and examples in .html format are included in the ./doc subdirectory.
    Help can also be accessed from Matlab via the `help' command.
#RUNNING
Start matlab from this directory or run startup.m to set the correct paths.

1. To run the scripts, first download the data and models from<br />
        <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/3DParallelBatching/data/> and <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/3DParallelBatching/models/>
       by typing `scons' in the ./data directory
2. The scripts that produced the examples shown in the documentation can be found in ./scripts.

#NOTES
In order to use the parallel version of 3D FWI, you need to open the parallel pool before running the script. If you want to use your local worker, 
just type the following command in the matlab command window:

`pool=parpool(b)` (n is the number of workers you want to use)

After the script finishes, just type the following command to close the parallel pool:

`delete(pool)` 


#SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
