%% 3D FWI framework example
%
% This software release includes a demonstration of the new 3D FWI framework developed here at SLIM. The goal 
% of this framework is to develop an FWI workflow that 
%  * reflects the mathematical framework of PDE constrained optimization as transparently as possible
%  * conforms to modern software engineering paradigms/designs
%  * is efficient, through the use of low-level operations implemented in C
%  * is easily parallelizable, through the use of Parallel Matlab (although other options are available)
%  * is easily comprehensible and extensible, so future features/optimizations can be seamlessly integrated into the codebase
%  * is tested and provably correct
%
% Author: Curt Da Silva
% Date  : July 2015


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/3DBasic/|.
%
% The code has been tested with _Matlab R2014b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/3DFreqModeling_ - 3D constant-density acoustic modeling
% * _solvers/Krylov  - iterative Krylov solvers
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.
% * _tools/solvers/QuasiNewton/minConf_mod/_ - nonlinear optimization

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/WaveformInversion/3DBasic/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% <example_edam1.html An example is shown here>.
