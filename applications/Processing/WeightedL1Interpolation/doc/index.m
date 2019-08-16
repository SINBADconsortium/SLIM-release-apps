%% Seismic Trace Interpolation using Weighted L1 Minimization 
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This software release includes a demonstration of the wavefield
% reconstruction framework in Matlab that has been developed at SLIM.
%
% Author: Hassan Mansour
%
% Date  : April, 2013


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Processing/WeightedL1Interpolation|.
%
% The code has been tested with _Matlab R2012b_.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _algorithms/AdaptiveSparseRecovery_ - weighted L1 minimization algorithm.
% * _operators/misc_ - Misc. operators such as source-receiver to midpoint-offset conversion.
% * _functions/misc_ - Misc. functions.
% * _solvers/SPGL1-SLIM_ - The spgl1-slim weighted L1 minimization solver.
% * _transforms/CurveLab-2.1.2-SLIM_ - curvelet transform functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/Processing/WeightedL1Interpolation/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed
% with the code.
% The scripts can be run in serial mode.

%% Scripts, examples and results
%
% Several scripts that reproduce the interpolation results in [1] can be found in the |exps| directory. 
%
% * |wL1_freq_SR|     - runs the reconstruction across frequency sclies in the source-receiver domain. 
% * |wL1_freq_MH|     - runs the reconstruction across frequency sclies in the midpoint-offset domain. 
% * |wL1_offset_TM|   - runs the reconstruction across offset sclies in the time-midpoint domain. 
%
% More detailed descriptions and some basic tests are described in <wL1_freq_MH.html wL1_freq_MH.m> and <wL1min.html wL1min.m>.
%
% <results.html The results are shown here>.

%% References
% 
% <https://slim.gatech.edu/Publications/Private/Submitted/Journal/mansour12iwr/mansour12iwr.pdf [1]> 
% H. Mansour, F. Herrmann, and O. Yilmaz. Improved wavefield reconstruction
% from randomized sampling via weighted one-norm minimization. submitted to
% Geophysics, 2013.

