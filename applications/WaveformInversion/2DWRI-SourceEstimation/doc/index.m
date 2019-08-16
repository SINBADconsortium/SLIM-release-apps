%% Waveform Reconstruction Inversion with Source Estimation
%
% This software release includes a demonstration of the  Wavefield Reconstruction
% Inversion technique with source estimation in Matlab that has been developed at SLIM.
%
% Author: Zhilong Fang
% Date  : March, 2018


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/2DWRI-SourceEstimation/|.
%
% The code has been tested with _Matlab R2015b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.
% * _tools/solvers/QuasiNewton/minConf_mod/_ - nonlinear optimization

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/WaveformInversion/2DWRI-SourceEstimation/|.
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
%
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% <example.html An example is shown here>.

%% References

%%
% <https://slim.gatech.edu/content/source-estimation-wavefield-reconstruction-inversion-0 [1]> Zhilong Fang, Rongrong Wang.
% Felix J. Herrmann, 2017, Source estimation for wavefield-reconstruction inversion.
