%% Waveform Reconstruction Imaging
%
% This software release includes a demonstration of the parallel Wavefield
% reconstruction Imaging algorithm in Matlab that has been developed at SLIM.
%
% Author: Bas Peters
% Date  : April, 2014


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Imaging/WRimaging.
%
% The code has been tested with _Matlab R2013a_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/2DFreqModeling_ - 2D constant-density acoustic modeling
% * _algorithms/WRI_ - Functions for the Wavefield Reconstruction Inversion algorithm
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/Imaging/WRimaging/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% <example_publish_Imaging.html An example is shown here>.

%% References

%%
% <http://dx.doi.org/10.1093/gji/ggt258 [1]> Tristan van Leeuwen, Felix J.
% Herrmann, Geophysical Journal International,2013. Mitigating local minima in full-waveform inversion by expanding the search space.

%%
% <https://slim.gatech.edu/content/penalty-method-pde-constrained-optimization
% [2]> Tristan van Leeuwen, Felix J. Herrmann. 2013. A penalty method for PDE-constrained optimization.

%%
% <https://slim.gatech.edu/content/examples-penalty-method [3]> Bas
% Peters, Felix J. Herrmann, Tristan van Leeuwen. EAGE, 2014. Wave-equation
% based inversion with the penalty method: adjoint-state versus wavefield-reconstruction inversion.
