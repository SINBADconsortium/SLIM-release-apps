%%Rank minimization based source-separation in time-jittered 2D ocean-bottom marine acquisition
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes a demonstration of the 2D time-jittered (or blended) marine 
% acquisition scheme, and a deblending framework based on rank minimization.
%
% Author: Rajiv Kumar (rakumar@eos.ubc.ca)
%
% Date: January, 2015


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Acquisition/2DTimeJitteredOBS-LR|.
%
% The code has been tested with _Matlab R2013a_ and do not requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _transforms/CurveLab-2.1.2-SLIM_ - CurveLab with SLIMs FFT Functions + other.
% * _tools/solvers/GenSPGL1_ - Generalized SPGL1.


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in
% this software release under |applications/Acquisition/2DTimeJitteredOBS-LR/examples/|. 
% Start matlab from |applications/Acquisition/2DTimeJitteredOBS-LR| to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed
% with the code.
% The scripts can be only run in serial mode.
% 

%% Functions
%
% The functions specific to this package can be found in the |misc_funcs| directory.


%% Examples and results
%
% <example.html Examples and results are shown here.>
% Scripts to reproduce the results can be found under
% |/applications/Acquisition/2DTimeJitteredOBS-LR/example|.


%% References
% <https://slim.gatech.edu/content/time-jittered-marine-acquisition-low-rank-vs-sparsity-0
% [1]> Rajiv Kumar, Haneet Wason, and F.J. Herrmann [2015] Time-jittered marine acquisition: low-rank v/s sparsity, EAGE Expanded Abstracts.
%
% <https://slim.gatech.edu/content/time-jittered-ocean-bottom-seismic-acquisition [2]> Haneet Wason, and F.J. Herrmann [2013] Time-jittered ocean bottom seismic acquisition, SEG Expanded Abstracts.
