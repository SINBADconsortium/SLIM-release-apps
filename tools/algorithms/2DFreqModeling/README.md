#2D Frequency-domain Full Waveform Inversion Environment in Matlab
This applications is available in the software release for both SINBAD
    consortium members and academic community.

##DESCRIPTION
This package contains Matlab functions for 2D frequency-domain modeling and waveform inversion.<br />

	F         - modeling operator
    DF        - Jacobian of modeling operator
    G         - Green's function for velocity linear profiles v0 + \alpha*z
    Helm2D    - Construct 2D Helmholtz operator as sparse matrix
    legendreQ - Legendre function, used for analytic Green's function
    oppDF     - pSPOT wrapper for DF
    test      - unit tests
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
Documentation for each of the functions can be accessed by typing `help <function>' in Matlab. 
    Examples are provided in applications/FWI 
##RUNNING
The functions can be called directly from Matlab. 
##NOTES
##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
