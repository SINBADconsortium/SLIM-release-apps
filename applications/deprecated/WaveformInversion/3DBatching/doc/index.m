%% 3D Frequency-domain FWI with batching
%
% This software release includes a demonstration of 3D frequency-domain
% FWI inversion using a row-projected Helmholtz solver described in 
%
% [1] T. van Leeuwen, D. Gordon, R. Gordon, and F.J. Herrmann, _Preconditioning the Helmholtz equation via row-projections_, in EAGE technical program, 2012
%
% [2] T. van Leeuwen, _Fourier analysis of the CGMN method for solving the Helmholtz equation_, arXiv:1210.2644
%
% [3] T. van Leeuwen and F.J. Herrmann, _3D Frequency-domain seismic inversion using a row-projected Helmholtz solver_, UBC Tech. report.
%
% Author: Tristan van Leeuwen
% Date  : March, 2013


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/3DBatching|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox for some of the functionality.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/3DFreqModeling_ - 3D variable-density acoustic modeling
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/3DBatching/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% Most scripts can be run in serial mode.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Functions
% The functions specific to this package can be found in the |mbin|
% directory.
% 
% * |misfit|      - Least-squares misfit function for 3D FWI including source-estimation
% * |lbfgs_batch| - L-BFGS with Wolfe linesearch and batching
% * |runFWI|      - driver for FWI
% * |runCARPCG|   - driver for rowprojected Helmholtz solver
% * |getoption|   - option parser

%% Examples and results
%
% Scripts to reproduce some of the examples found in [3] can be found in the |scripts| directory. 
% <results.html The results are shown here>.


