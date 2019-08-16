%% Large scale, parallel low rank matrix completion for seismic data interpolation
%
% This applications is available only in the software release for members of SINBAD consortium.
%
% This software implements the SPGLR algorithm using block parallel matrix-matrix 
% multiplication, which allows it to be used on large scale problems.
%
% Author: Curt Da Silva (curtd@math.ubc.ca)
%
% Date: February, 2015

%% Downloading & Dependencies
% The synthetic examples code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |applications/Processing/LargeScaleLRI|.
%
% The code has been tested with _Matlab R2013a_ and requires the Parallel Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part of the 
% <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT-SLIM      - object oriented framework for matrix-free linear algebra.
% * _utilities/functions/misc - miscellaneous functions, for plotting, experiment organization, etc.
% * _tools/solvers/SPGLR_PAR  - parallel SPGLR package


%% Running & Parallelism
% All of the examples and results are produced by the scripts found in this software release under |/applications/Processing/LargeScaleLRI/examples|.
% Start matlab from |/applications/Processing/LargeScaleLRI| to add the appropriate paths.
%
% To run the scripts, follow the instructions in the README.md file enclosed with the code

%% Functions
% The SPGLR code consists of the spgLR.m function |_tools/solvers/SPGLR_PAR|, which performs the parameter-cooling method by solving the Lasso subproblems, and spgLRobj.m, which performs the distributed computation of the objective and its gradient. The accompanying documentation is in the README.md file.
% 

%% Examples and results
% The examples of large scale missing trace interpolation using these methods can be found in 
% |applications/Processing/LargeScaleLRI/examples|
% 
% Results are of missing-receiver interpolation is shown in <spgLR_bgdata_view.html spgLR_bgdata_view.m>

%% References
% <https://slim.gatech.edu/content/fast-methods-denoising-matrix-completion-formulations-applications-robust-seismic-data-inter [1]> A. Aravkin, et al. "Fast methods for denoising matrix completion formulations, with applications to robust seismic data interpolation"

%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.
