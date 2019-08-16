A much nicer-formatted online README along with a short execution walkthrough
is available here:<br />
<https://slim.gatech.edu/SoftwareDemos/applications/WavefieldSeparation/RobustEPSI/><br />
It is a good idea to refer to the above URL first and only refer to this
'README.txt' and 'README_examples.txt' files for detailed technical information.

#Robust Estimation of Primaries with Sparse Inversion (via L1 minimization) for Pre-Stack Seismic Lines


##DESCRIPTION

This package is an application for simultaneous multiple removal and
    estimation of the source signature. It is based on the Estimation of
    Primaries by Sparse Inversion approach but reformulated via block-
    coordinate descent and convexification to L1 minimization. We believe this
    is therefore more robust to artifacts, and should converge quicker and
    more reliably that the original formulation.

##COPYRIGHT

You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.

##PREREQUISITES

1. MATLAB r2009b (7.9) or above installed on your system<br />
This software package is almost entirely written in Matlab. Make sure 
        you can successfully invoke the Matlab prompt and display figures.<br />
        To use the standalone program, make sure the command-line interface to 
        Matlab can be invoked from your shell with the command 'matlab 
        -nodesktop'. Usage of the parallel version requires the Parallel 
        Computing Toolbox (v4.2+).
    
2. SPOT (provided with this release)<br />
The SPOT operator framework is used to help keep our code modular and readable.
        
3. (optional) CurveLab 2.1.2 w/Matlab interface compiled on your system<br />
(provided with this release)<br />
        In order to turn on sparsity transform features in our software demo 
        and program, CurveLab and its Matlab interface must be compiled and 
        installed on your system. Make sure Matlab can find the compiled 
        CurveLab mex files.<br />
        * Run the included 'test_sparsityTransform.m' script to test 
        whether all relevant sparsity transform programs are installed 
        correctly.
        
4. Python 2.3 or above installed on your system (http://www.python.org)<br />
        SCons 1.2.0 or above installed on your system (<http://www.scons.org>)<br />
        cURL 7 or above installed on your system (<http://curl.haxx.se/>)<br />
        These programs are used in the demo to construct scripts and fetch 
        datafiles from an FTP server.
        
NOTE: Most prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
        
##INSTALLATION NOTES
Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.<br />
Before using this program you should make sure your environment is setup
    correctly. Make sure you have the 'SLIM_ROOT' environment variable set to
    the directory containing the SLIM software release. We recommend also
    sourcing the 'environment.sh/csh' file in your shell login script, like
    so (for BASH):
            
            CURR_DIR=`pwd`
            cd /path/to/slim
            . environment.sh
            cd $CURR_DIR
            
Run the included 'test_SPOT.m' script from your Matlab prompt to test
    whether SPOT has been installed properly.

Run the included 'test_sparsityTransforms.m' script to test whether all
    relevant sparsity transform programs are installed correctly.
	
In order to use the standalone programs, the 'bin' directory should be 
    added to your shell's search path.
    
##RUNNING EXAMPLES

1. Start matlab from this directory or run startup.m to set the correct
      paths.
2. Consult the 'README_examples.txt' file for more information.

##STANDALONE PROGRAM
The Stabilized EPSI program can be run as standalone command-line shell 
    programs instead of a demo. The executables, located the 'bin' directory, 
    are executable Python scripts. Consult the 'robust_epsi_README.txt' file 
    in there for more information on these programs.

##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.


##REFERENCES
The original Estimation of Primaries by Sparse Inversion formulation:

[1] G.J. van Groenestijn and D.J. Verschuur. Estimating primaries by 
        sparse inversion and application to near-offset data reconstruction. 
        Geophysics (2009) vol. 74 (3) pp. A23-A28
        
Estimation of Primaries with L1 Inversion:
    
[2] T.T.Y. Lin and F.J. Herrmann.  Unified compressive sensing framework 
        for simultaneous acquisition with primary estimation. SEG Technical 
        Program Expanded Abstracts (2009) vol. 28 pp. 3113-3117
    
[3] T.T.Y. Lin and F.J. Herrmann. Stabalized Estimation of Primaries by 
        Sparse Inversion. EAGE Technical Program Expanded Abstracts (2010).
        
[4] T.T.Y. Lin and F.J. Herrmann. Estimating Primaries by Sparse Inversion
        in a Curvelet-like Representation Domain. EAGE Technical Program 
        Expanded Abstracts (2011).


