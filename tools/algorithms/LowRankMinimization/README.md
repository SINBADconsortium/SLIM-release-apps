#Seismic data regularization, interpolation and denoising using SVD-free low-rank matrix factorization.

##DESCRIPTION
This package contains Matlab functions for seismic data interpolation and denoising.
    runinterp         - read all the user sepcified input data and parameters
    timetofreq        - Convert input time domain data into frequency domain to
                        perform interpolation and/ or denoising on monochromatic
                        frequency slices. We only select the positive part of spectra.
    freqtotime        - Convert interpolated and/or denoised output data into time domain.
    syndata           - Make data set for test case attached with the package.
    opMH              - Perform source-receiver to midpoint-offset transformation on each
                        monochromatic frequency slices.
    LowRank_2D        - Distribute the input data and call main function one each
                        monochromatic frequency slices to perform interpolation
                        and/ or denoising.
    Interp_Denoise    - Main function to perform the interpolation and /or denoising
                        on each monochromatic frequency slices.
    Reg_Interp        - Main function to perform the regularization, interpolation and /or denoising
                        on each monochromatic frequency slices.

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
The functions can be called directly from Matlab.

##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##REFERENCES
 [1] A.Y. Aravkin, R. Kumar, H. Mansour, B. Recht, F. J. Herrmann, 2014. Fast methods for denoising matrix completion formulations, with application to robust seismic data interpolation.

 [2] R. Kumar, O. Lopez, E. Esser, F. J. Herrmann, 2014. Matrix completion on unstructured grids : 2-D seismic data regularization and interpolation, submitted to SEG.

 [3] R. Kumar, A.Y. Aravkin, H. Mansour, B. Recht, F. J. Herrmann, 2013. Seismic data interpolation and denoising using SVD-free low-rank matrix factorization, EAGE.

