#Seismic data regularization, interpolation and denoising using SVD-free low-rank matrix factorization.

##DESCRIPTION
This package illustrates the main features of the regularization, missing-trace 
    interpolation and denoising. This package includes 1) the algorithm 
    2) synthetic examples. The algorithm is in directory "function" 
    and the examples are in "examples",both under current folder.

##ON-LINE DOCUMENTATION
<https://slim.gatech.edu/SoftwareDemos/applications/Processing/LowRankInterpolationAndDenoising/>

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
Start matlab from this directory or run startup.m to set the correct paths. The scripts 
     that produced the examples can be found in ./examples/GulfofSuez. The scripts can be run on
     either serial or parallel version of matlab. Set "options.parallel" in each example file
     according to the requirements 

1. To run the "GofS_Interp_and_denoise.m" and "GofS_Interp.m" scripts in example folder, download the data from
    	 <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Common/GulfOfSuez>
        by typing `scons' in the ./data/Interpolation directory

2. To run the "GofS_regularize_and_interp.m" scripts in example folder, download the data from
    	 <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Processing/LowRankInterpolationAndDenoising/Regularization>
        by typing `scons' in the ./data/Regularization directory

2. You can find three sets of examples under "./examples/GulfofSuez" : First example
        is "GofS_Interp", where we demonstrate the missing-trace interpolation,
        second example is "GofS_Interp_and_denoise", where we demostrate
        missing-trace interpolation and denoising and third example is "GofS_regularize_and_interp.m"
        where we domenostarte regularization.  

3. The output results from example section will be stored in "./results"
        directory. The algorithm will generate the ".dat" file for each of the freqeuency
        slices indexed by frequency number, which will contain the SPGL1 log
        information.

4. Running your own data<br />
    4.1 Note: To run the algorithm on your own data set requires 
          equal number of sources and receivers in input data set.

5. The current version of algorithm works only with parallel matlab.

##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


##REFERENCES
 [1] A.Y. Aravkin, R. Kumar, H. Mansour, B. Recht, F. J. Herrmann, 2014. Fast methods for denoising matrix completion formulations, with application to robust seismic data interpolation.

 [2] R. Kumar, O. Lopez, E. Esser, F. J. Herrmann, 2014. Matrix completion on unstructured grids : 2-D seismic data regularization and interpolation, submitted to SEG.

 [3] R. Kumar, A.Y. Aravkin, H. Mansour, B. Recht, F. J. Herrmann, 2013. Seismic data interpolation and denoising using SVD-free low-rank matrix factorization, EAGE.
