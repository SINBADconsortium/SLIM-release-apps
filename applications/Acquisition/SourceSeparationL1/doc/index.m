%% Source separation for towed-streamer marine data via sparsity promotion
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This package contains a MATLAB implementation of a 2-D over/under blended marine acquisition scheme, and a deblending (or source separation) algorithm based on sparsity-promoting inversion in the curvelet domain using L1 minimization.
% 
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: July, 2016


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under |applications/Acquisition/SourceSeparationL1|.
%
% The code has been tested with _Matlab R2015b_.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _functions/misc_ - miscellaneous functions.
% * _solvers/SPGL1-SLIM_ - SLIM version of SPGL1 solver.
%
% The CurveLab package is downloaded from the |tools| part of the SLIM-release-comp repository. 
% * _transforms/CurveLab-2.1.2-SLIM_ - curvelet transform functions.


%% Functions
%
% Some functions specific to this package can be found in the |applications/Acquisition/SourceSeparationL1/misc_funcs| directory.


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in this software release under |applications/Acquisition/SourceSeparationL1|. Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file under |applications/Acquisition/SourceSeparationL1/examples|. The scripts are run in serial mode but can also be run in parallel with slight modifications to the code. Please see the README file.
%
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.


%% Examples and results
%
% <example.html Examples and results are shown here.> Scripts to reproduce the results can be found under |applications/Acquisition/SourceSeparationL1/examples|. 


%% References
%
% <https://slim.gatech.edu/content/source-separation-simultaneous-towed-streamer-marine-acquisition-–-compressed-sensing-approa [1]> Rajiv Kumar, Haneet Wason, and Felix J. Herrmann [2015]. Source separation for simultaneous towed-streamer marine acquisition --- a compressed sensing approach, Geophysics, vol. 80, pp. WD73-WD88.


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

