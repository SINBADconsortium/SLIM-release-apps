%% Modified gauss-newton full-waveform inversion
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes an parallel framework in Matlab for modified 
% gauss-newton (GN) full-waveform inversion [1,2], which based on the ideas from 
% compressive-sensing and stochastic optimization, where the GN updates 
% are computed from random subsets of the data via curvelet-domain 
% sparsity-promotion. There are two different subset sampling strategies 
% are considered in this package: randomized source encoding, and drawing 
% sequential shots firing at random source locations from marine data 
% with missing near and far offsets. In both cases, we obtain excellent 
% inversion results compared to conventional methods at reduced 
% computational costs. There is also a small example based on Camembert
% example which can allow users to test the algorithm in a short time 
%
% Author: Xiang Li
% Date  : March, 2012

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware/ SLIM sofware release> under
% |/applications/WaveformInversion/2DModGaussNewton|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware/ SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/2DFreqModeling_ - 2D constant-density acoustic modeling
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.


%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/2DModGaussNewton/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% Scripts to reproduce the famous Camembert example [7],  as well as
% results from sevaral papers can be foundin the |scripts| directory. <examples.html
% The results are shown here>.
%

%% References
%
% <https://slim.gatech.edu/node/6390 [1]> Felix J. Herrmann, Xiang Li, Aleksandr Y. Aravkin, and Tristan van Leeuwen, A modified, sparsity promoting, Gauss-Newton algorithm for seismic waveform inversion, in Proc. SPIE, 2011, vol. 2011.
%
% <https://slim.gatech.edu/node/6621 [2]> Xiang Li, Aleksandr Y. Aravkin, Tristan van Leeuwen, and Felix J. Herrmann, Fast randomized full-waveform inversion with compressive sensing. 2011. Geophysics, accepted.
%
% <http://dx.doi.org/10.1190/1.1442188  [3]> O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)
%
%% Acknowledgements
% The synthetic Compass model was provided by the BG-GROUP, see also the
% <BG_DISCLAIMER.txt disclaimer>.
