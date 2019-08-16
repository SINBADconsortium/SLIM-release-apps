%% 2D ocean-bottom marine acquisition via jittered sampling
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes a demonstration of the 2D time-jittered (or blended) marine 
% acquisition scheme, and a deblending framework based on sparse inversion via L1 minimization.
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: April, 2013


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Acquisition/2DTimeJitteredOBS|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _transforms/CurveLab-2.1.2-SLIM_ - CurveLab with SLIMs FFT Functions + other.
% * _solvers/SPGL1-SLIM_ - SLIM version of SPGl1


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in
% this software release under |applications/Acquisition/2DTimeJitteredOBS/examples/|. 
% Start matlab from |applications/Acquisition/2DTimeJitteredOBS| to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed
% with the code.
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.


%% Functions
%
% The functions specific to this package can be found in the |acq_funcs| directory.


%% Examples and results
%
% <examples.html Examples and results are shown here.>
% Scripts to reproduce the results can be found under
% |/applications/Acquisition/2DTimeJitteredOBS/examples|.


%% References
%
% <https://slim.gatech.edu/content/time-jittered-ocean-bottom-seismic-acquisition [1]> Haneet Wason, and F.J. Herrmann [2013] Time-jittered ocean bottom seismic acquisition, SEG Expanded Abstracts.
%
% <https://slim.gatech.edu/content/ocean-bottom-seismic-acquisition-jittered-sampling-0 [2]> Haneet Wason, and F.J. Herrmann [2013] Ocean bottom seismic acquisition via jittered sampling, EAGE Expanded Abstracts.
%
% <https://slim.gatech.edu/content/randomized-marine-acquisition-compressive-sampling-matrices [3]> Hassan Mansour, Haneet Wason, Tim T.Y. Lin, and Felix J. Herrmann [2012] Randomized marine acquisition with compressive sampling matrices, Geophysical Prospecting, vol. 60, 648-662.

