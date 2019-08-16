%% Basic Waveform inversion
%
% This software release includes a demonstration of the parallel waveform
% inversion framework in Matlab that has been developed at SLIM.
%
% Author: Tristan van Leeuwen
% Date  : May, 2012


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/2DBasic|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/2DFreqModeling_ - 2D constant-density acoustic modeling
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/2DBasic/|. 
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
% * |Jls|        - least-squares residual and gradient for given model and observed data.
% * |mylbfgs|    - standard L-BFGS with Wolfe linesearch
%
% <mfiles.html More detailed descriptions and some basic tests are described here>.

%% Examples and results
%
% A script to reproduce the famous Camembert example [1] can be found in the |scripts| directory. 
% <results.html The results are shown here>.

%% References
% 
% <http://dx.doi.org/10.1190/1.1442188 [1]> O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear
% inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)

