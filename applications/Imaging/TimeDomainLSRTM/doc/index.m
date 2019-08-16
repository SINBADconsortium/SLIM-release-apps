%% Time domain LSRTM with sparsity promotion
% This applications is available only in the software release for
% members of SINBAD consortium.
%
% This software provides an algorithm to perform least-squares RTM with
% sparsity promotion using the linearized Bregman method. By subsampling
% the sources in each iteration, the overall number of PDE solves lies in
% the range of a regular RTM image, making this method feasible for
% large-scale problems.
%
% Author: Philipp Witte (pwitte@eos.ubc.ca)
% Date: February 2016

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |/applications/Imaging/TimeDomainLSRTM|
%
% The code has been tested with _Matlab R2015a_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/TimeModeling_ - 2/3D acoustic modeling
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.
%
% If you want to use your own modules to do modelling or multiple
% prediction, please contact the author.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/Imaging/TimeDomainLSRTM/examples/|. 
% Start matlab from |/applications/Imaging/TimeDomainLSRTM| to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file
% enclosed inside the folder for each set of examples.
% 
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
% <example.html Examples and results are shown here.>
% Scripts to reproduce the results can be found under
% |/applications/Imaging/TimeDomainLSRTM/examples|.

