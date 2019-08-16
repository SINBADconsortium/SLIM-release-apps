%% Fast Robust Waveform inversion
%
% This software release includes a demonstration of the parallel waveform
% inversion framework in Matlab that has been developed at SLIM. In particular, it
% includes implementations of waveform inversion using robust penalties and
% source estimation [1,2] in conjunction with an optimization strategy that
% works with only a few sources at each iteration [3-6]. Updated to conform to the new
% waveform inversion software framework [8].
%
% Author: Tristan van Leeuwen, Aleksandr Aravkin
% Date  : March, 2012
%
% Updated : Curt Da Silva
% Date : May, 2016


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/2DRobustBatching|.
%
% The code has been tested with _Matlab R2015b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/2DFreqModeling_ - 2D constant-density acoustic modeling
% * _algorithms/CommonFreqModeling - baseline acoustic modeling components
% * _solvers/QuasiNewton - quasi newton optimization algorithms
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/2DRobustBatching/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Functions
% The functions specific to this package can be found in the |mbin|
% directory.
% 
% * |twonorms|   - the two-norm and its derivative. 
% * |hubers|     - the Huber misfit and its derivative. 
% * |hybrid|     - the hybrid misfit and its derivative.
% * |students|   - the Students t misfit and its derivative.
% * |mylbfgs|    - standard L-BFGS with Wolfe linesearch
% * |minfunc_batch| - stochastic optimization with bound constraints
%
% <mfiles.html More detailed descriptions and some basic tests are described here>.

%% Examples and results
%
% Scripts to reproduce the famous Camembert example [7],  as well as
% results from several papers can be found in the |scripts| directory. <results_2drobustbatching.html
% The results are shown here>.

%% References
% 
% <https://slim.gatech.edu/node/6610 [1]> A.Y. Aravkin, T. van Leeuwen and F.J. Herrmann - Source estimation
% for frequency-domain FWI with robust penalties, EAGE Expanded abstracts
% 2012.
%
% <http://dx.doi.org/10.1190/1.3627747 [2]> A.Y. Aravkin, T. van Leeuwen, and F.J. Herrmann, 2011. Robust full-waveform inversion using the Student's t-distribution. SEG Technical Program Expanded Abstracts 30, 2669-2673.
%
% <http://arxiv.org/abs/1104.2373 [3]> M.P. Friedlander, M.Schmidt, 2011. Hybrid Deterministic-Stochastic Methods for Data Fitting. arXiv:1104.2373.
%
% <http://www.earthdoc.org/detail.php?pubid=50341 [4]> T. van Leeuwen, F.J. Herrmann, M. Schmidt, M.P. Friedlander, 2011. A hybrid stochastic-deterministic optimization method for waveform inversion. EAGE Expanded Abstracts.
%
% <https://slim.gatech.edu/node/6620 [5]> T. van Leeuwen and F.J. Herrmann - Fast Waveform inversion without
% source-encoding, Geophysical Prospecting, submitted
%
% <https://slim.gatech.edu/node/6622 [6]> A.Y. Aravkin, M.P. Friedlander, F.J. Herrmann, and T. van Leeuwen, 2011. Robust inversion, dimensionality reduction, and randomized sampling. Mathematical Programming, accepted.
%
% <http://dx.doi.org/10.1190/1.1442188 [7]> O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear
% inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)
%
% <https://slim.gatech.edu/content/unified-2d3d-software-environment-large-scale-time-harmonic-full-waveform-inversion [8]> C. Da Silva, F.J. Herrmann - A unified 2D/3D software environment for large scale time-harmonic full waveform inversion.

%% Acknowledgements
% The synthetic Compass model was provided by the BG-GROUP, see also the
% <BG_DISCLAIMER.txt disclaimer>.
