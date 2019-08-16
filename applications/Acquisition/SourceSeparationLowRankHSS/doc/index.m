%% Source separation via SVD-free rank minimization in the hierarchical semi-separable representation
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes a demonstration of a 2-D over/under blended marine acquisition scenario, and a source separation algorithm based on a SVD-free rank minimization scheme.
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: June, 2014; July, 2016 - parallel version of the code made available.

%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under |/applications/Acquisition/SourceSeparationLowRankHSS|.
%
% The code has been tested with _Matlab R2015b_. The updated code requires the Parallel Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT-SLIM_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _solvers/GenSPGL1_ - Generalized SPGL1.


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in this software release under |applications/Acquisition/SourceSeparationLowRankHSS/examples/|. Start matlab from |applications/Acquisition/SourceSeparationLowRankHSS| to add the appropriate paths. 
% 
% To run the scripts follow the instructions in the README file enclosed with the code. The scripts are run in parallel mode.


%% Examples and results
%
% <example.html Examples and results are shown here.> 
% Scripts to reproduce the results can be found under |/applications/Acquisition/SourceSeparationLowRankHSS/examples|.


%% References
%
% <https://slim.gatech.edu/content/source-separation-simultaneous-towed-streamer-marine-acquisition-–-compressed-sensing-approa [1]> Rajiv Kumar, Haneet Wason, and Felix J. Herrmann [2015]. Source separation for simultaneous towed-streamer marine acquisition --- a compressed sensing approach, Geophysics, vol. 80, pp. WD73-WD88.
%
% <https://slim.gatech.edu/content/source-separation-svd-free-rank-minimization-hierarchical-semi-separable-representation [2]> Haneet Wason, Rajiv Kumar, Aleksandr Y. Aravkin, and Felix J. Herrmann [2014]. Source separation via SVD-free rank minimization in the hierarchical semi-separable representation, presented at the SEG Conference.
%
% <https://slim.gatech.edu/content/fast-methods-denoising-matrix-completion-formulations-applications-robust-seismic-data-inter [3]> Aleksandr Y. Aravkin, Rajiv Kumar, Hassan Mansour, Ben Recht, and Felix J. Herrmann [2014]. Fast methods for denoising matrix completion formulations, with applications to robust seismic data interpolation, SIAM Journal on Scientific Computing (SISC), vol. 36, pp. S237-S266.
%
% <https://slim.gatech.edu/content/reconstruction-seismic-wavefields-low-rank-matrix-factorization-hierarchical-separable-matri [4]> Rajiv Kumar, Hassan Mansour, Aleksandr Y. Aravkin, and Felix J. Herrmann [2013]. Reconstruction of seismic wavefields via low-rank matrix factorization in the hierarchical-separable matrix representation, SEG Technical Program Expanded Abstracts, p. 3628-3633.


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

