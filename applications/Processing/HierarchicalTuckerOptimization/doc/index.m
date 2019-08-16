%% Missing trace interpolation of 3D seismic data using the
%% Hierarchical Tucker tensor format
%
% This applications is available only in the software release for members of SINBAD consortium.
%
% This software provides an algorithm for missing trace/receiver
% interpolation of 3D seismic data using the latent Hierarchical
% Tucker tensor format. The algorithm operates on a single
% frequency slice, although an extension to an arbitrary number of
% frequency slices is straightforward.  
%
% Author: Curt Da Silva (curtd@math.ubc.ca)
%
% Date: March, 2014

%% Downloading & Dependencies
% The synthetic examples code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |applications/Processing/HierarchicalTuckerOptimization|.
%
% The code has been tested with _Matlab R2013a_ and supports, but does not require, the Parallel Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part of the 
% <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT-SLIM - object oriented framework for matrix-free linear algebra.
% * _tools/solvers/HTOpt - Hierarchical Tucker tensor algorithms

%% Running & Parallelism
% All of the examples and results are produced by the scripts found in this software release under |/applications/Processing/HierarchicalTuckerOptimization/examples|.
% Start matlab from |/applications/Processing/HierarchicalTuckerOptimization| to add the appropriate paths.
%
% To run the scripts, follow the instructions in the README.md file enclosed with the code

%% Functions
% The missing-trace tensor interpolation using Hierarchical Tucker
% tensors can be found in |_tools/algorithms/HTOpt|, with
% accompanying documentation in the README.md file.
% 

%% Examples and results
% Several examples (serial and parallel versions) of missing trace
% interpolation using these methods can be found in 
% |applications/Processing/HierarchicalTuckerOptimization/examples|
% 
% Results are of missing-receiver interpolation is shown in <interp4Dview.html interp4Dview.m>

%% References
% <https://slim.gatech.edu/Publications/Public/TechReport/2014/dasilva2014htuck/dasilva2014htuck.pdf [1]> C. Da Silva and F. J. Herrmann, 2014. Optimization on the Hierarchical Tucker manifold - applications to tensor completion
%
%
% <https://slim.gatech.edu/Publications/Public/Conferences/SAMPTA/2013/dasilva2013SAMPTAhtuck/dasilva2013SAMPTAhtuck.pdf [2]> C. Da Silva and F. J. Herrmann, 2013. Hierarchical Tucker Tensor Optimization - Applications to Tensor Completion

%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.
