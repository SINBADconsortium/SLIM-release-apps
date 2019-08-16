%% Parallel 3D frequency domain full waveform inversion
%
% This software release includes a demonstration of 3D frequency-domain
% FWI inversion using a stencil-based Helmholtz matrix-vector multiply, Multi-level GMRES preconditioner, automatic parallelization over source/frequency, and a stochastic LBFGS scheme, described in
% 
% [1] "A unified 2D/3D software environment for large scale time-harmonic full waveform inversion", Curt Da Silva and Felix Herrmann, SEG 2016 (Submitted).
%
% Author: Curt Da Silva
% Date  : April, 2016


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/3DParallelBatching|.
%
% The code has been tested with _Matlab R2015b_ and requires the Parallel
% Computing Toolbox for some of the functionality.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/3DFreqModeling_ - 3D time-harmonic acoustic modeling
% * _algorithms/CommonFreqModeling_ - general tools for 2D/3D time-harmonic FWI
% * _solvers/Krylov - linear system methods + options
% * _solvers/Multigrid - Multigrid preconditioner for the Helmholtz equation
% * _solvers/QuasiNewton - Quasi Newton optimization methods
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/3DParallelBatching/|. 
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
% * |dload_fslices|   - Load frequency slices in a distributed fashion
% * |fwi_exp##|       - Parameter files corresponding to different FWI experiments
% * |fwi_exp|         - Main driver code for running inversion with different experiment settings
% * |getoption|       - Option parsing
% * |load_geometry|   - Loads the acquisition geometry corresponding to a particular model
% * |load_model|      - Loads the velocity associated to a particular model
% * |minfunc_frugal|  - Stochastic optimization with various sub-problem optimization kernels
% * |plot_fwi_exp|    - Results plotting
% * |rand_set|        - Get a binned random vector
% * |vel_plot|        - Velocity plotting


%% Examples and results
%
% Scripts to generate results of examples can be found in the |scripts| directory. 
% <results.html The results are shown here>.


