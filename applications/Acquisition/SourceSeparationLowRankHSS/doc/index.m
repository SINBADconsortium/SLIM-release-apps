%% Source separation via SVD-free rank minimization in the hierarchical semi-separable representation
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes a demonstration of a 2-D over/under, blended marine acquisition scenario, and a source separation algorithm based on a SVD-free rank minimization scheme.
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: June, 2014


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under |/applications/Acquisition/SourceSeparationLowRankHSS|.
%
% The code has been tested with _Matlab R2013a_ and supports, but does not require, the Parallel Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT-SLIM_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _solvers/GenSPGL1_ - Generalized SPGL1


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in this software release under |applications/Acquisition/SourceSeparationLowRankHSS/examples/|. Start matlab from |applications/Acquisition/SourceSeparationLowRankHSS| to add the appropriate paths. 
% 
% To run the scripts follow the instructions in the README file enclosed with the code. The scripts can be run in serial mode.


%% Examples and results
%
% <example.html Examples and results are shown here.> 
% Scripts to reproduce the results can be found under |/applications/Acquisition/SourceSeparationLowRankHSS/examples|.


%% References
%
% <https://slim.gatech.edu/content/source-separation-svd-free-rank-minimization-hierarchical-semi-separable-representation [1]> Haneet Wason, Rajiv Kumar, Aleksandr Y. Aravkin, and F.J. Herrmann [2014]. Source separation via SVD-free rank minimization in the hierarchical semi-separable representation, to be presented at the SEG Conference.
%
% <https://slim.gatech.edu/content/fast-methods-denoising-matrix-completion-formulations-application-robust-seismic-data-interp [2]> Aleksandr Y. Aravkin, Rajiv Kumar, Hassan Mansour, Ben Recht, and F.J. Herrmann [2013]. Fast methods for denoising matrix completion formulations, with application to robust seismic data interpolation, submitted to the SIAM Journal on Scientific Computing (SISC).
%
% <https://slim.gatech.edu/content/reconstruction-seismic-wavefields-low-rank-matrix-factorization-hierarchical-separable-matri [3]> Rajiv Kumar, Hassan Mansour, Aleksandr Y. Aravkin, and Felix J. Herrmann [2013]. Reconstruction of seismic wavefields via low-rank matrix factorization in the hierarchical-separable matrix representation, SEG Technical Program Expanded Abstracts, p. 3628-3633.


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

