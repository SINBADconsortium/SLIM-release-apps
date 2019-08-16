%% Seismic data regularization, missing-trace interpolation and denoising using SVD-free low-rank matrix factorization.
%
%
% This applications is available only in the software release for members of SINBAD consortium.
%
%
% This software provides an algorithm for seismic data regularization,
% missing-trace interpolation and denoising (using Generalized SPGl1 as solver). The algorithm
% solves the system in parallel over frequencies. 
% The regularization, missing-trace interpolation and denoising is done using 
% robust-rank regularized formulation. We illustrate the advantages of 
% the new approach using a seismic line from Gulf of Suez.
%
% Author: Rajiv Kumar (rakumar@eos.ubc.ca)
%
% Date: April,2014

%% Downloading & Dependencies
% The synthetic examples code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |applications/Processing/LowRankInterpolationAndDenoising|.
%
% The code has been tested with _Matlab R2013a_ and require the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT-SLIM_ - object oriented framework for matrix-free linear algebra.
% * _tools/algorithms/LowRankMinimization_ - Matrix factorization based
%                                            low-rank optimization algorithm.
% * _tools/solvers/GenSPGL1_ - Generalized SPGL1.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/Processing/LowRankInterpolationAndDenoising/examples/|. 
% Start matlab from |/applications/Processing/LowRankInterpolationAndDenoising| to add the appropriate paths.
%
% 
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
%% Functions
% The regularization, missing-trace interpolation and denoising code can be found in 
% |_tools/algorithms/LowRankMinimization_|. The main components are listed below 
%
% _algorithms/LowRankMinimization_
%
% * |runinterp|  - Read the input parameter and data 
% * |opMH|       - Transform monochromatic frequency slice from
%                  the source-receiver domain to the midpoint-offset domain.
% * |LowRank_2D| - Perform regularization, missing-trace interpolation and denoising in
%                  the midpoint-offset domain.
%

%% Examples and results
% An examples of regulation, missing-trace interpolation and denoising can be found in 
% |applications/Processing/LowRankInterpolationAndDenoising|
%
% Results of regularization, missing-trace interpolation and denoising is shown in <GofS_Interp.html GofS_Interp.m>.

%% References
% <https://slim.gatech.edu/Publications/Private/Conferences/SEG/2014/kumar2014SEGmcu/kumar2014SEGmcu.html [1]> R. Kumar, O. Lopez, E. Esser, F. J. Herrmann, 2014. Matrix completion on unstructured grids : 2-D seismic data
% regularization and interpolation, submitted to SEG. 
%
% <http://arxiv.org/abs/1302.4886 [2]> A.Y. Aravkin, R. Kumar, H. Mansour,
% B. Recht, F. J. Herrmann, 2013. Fast methods for denoising matrix completion formulations, with application to robust seismic data interpolation..
%
% <https://slim.gatech.edu/Publications/Public/Conferences/EAGE/2013/kumar2013EAGEsind/kumar2013EAGEsind.pdf [3]> R. Kumar, A.Y. Aravkin, H. Mansour,
% B. Recht, F. J. Herrmann, 2013. Seismic data interpolation and denoising using SVD-free low-rank matrix factorization, EAGE.

%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.
