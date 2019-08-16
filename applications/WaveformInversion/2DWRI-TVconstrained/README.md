# Total Variation Regularized Wavefield Reconstruction Inversion

##  DESCRIPTION

This MATLAB code implements a scaled gradient projection method to minimize the wavefield reconstruction inversion (WRI) objective subject to total variation and spatially varying bound constraints.  It extends an earlier implementation of the penalty method for full wave form inversion by Tristan van Leeuwen and complements the 2DWRI implementation at <a href = "https://github.com/SINBADconsortium/SLIM-release-apps/tree/master/applications/WaveformInversion/2DWRI">https://github.com/SINBADconsortium/SLIM-release-apps/tree/master/applications/WaveformInversion/2DWRI</a>.  

##  ON-LINE DOCUMENTATION

Algorithm details can be found in the technical report at <a href="https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion">https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion</a> and the SINBAD presentation at <a href="https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion-0">https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion-0</a>  A summary of the included MATLAB code can be found at <a href="https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI-TVconstrained">https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI-TVconstrained</a>.

##  COPYRIGHT
 
You may use this code only under the conditions and terms of the license contained in the files LICENSE or COPYING provided with this source code. If you do not agree to these terms you may not use this software.

##  PREREQUISITES

All prerequisites, except for MATLAB, are provided in the software release repositories and should be installed as necessary before using any of SINBAD's software.  Although the current code requires only MATLAB, future versions will require functions defined in the software release repositories.

##  INSTALLATION

###  Software in SLIM-release-comp repository
 
This application does not require any installation of extra software from SLIM-release-comp repository.

##  RUNNING

###  Downloading data

The part of the BP 2004 velocity benchmark data set used as the true velocity model for one of the included examples can be downloaded from <a href="ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/2DWRI-TVconstrained">ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/2DWRI-TVconstrained</a>.

###  Running applications/demos

Examples can be run by typing main.m.  The data, model parameters and algorithm parameters are all defined in setup_problem.m.  See <a href="https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI-TVconstrained">https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/2DWRI-TVconstrained</a> for more details.

####  Hardware requirements

This code was developed using MATLAB R2013a.  

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES

[1] http://dx.doi.org/10.1093/gji/ggt258 Tristan van Leeuwen, Felix J. Herrmann, Geophysical Journal International,2013. Mitigating local minima in full-waveform inversion by expanding the search space.

[2] https://slim.gatech.edu/content/penalty-method-pde-constrained-optimization Tristan van Leeuwen, Felix J. Herrmann. 2013. A penalty method for PDE-constrained optimization.

[3] https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion Ernie Esser, Tristan van Leeuwen, Aleksandr Aravkin and Felix Herrmann.  Technical Report, 2014.  A scaled gradient projection method for total variation regularized full waveform inversion

[4] https://slim.gatech.edu/content/scaled-gradient-projection-method-total-variation-regularized-full-waveform-inversion-0 Ernie Esser. SINBAD Meeting, Spring 2014.  A scaled gradient projection method for total variation regularized full waveform inversion

[5] https://slim.gatech.edu/content/examples-penalty-method Bas Peters, Felix J. Herrmann, Tristan van Leeuwen. EAGE, 2014. Wave-equation based inversion with the penalty method: adjoint-state versus wavefield-reconstruction inversion.
