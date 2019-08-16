%% Robust Estimation of Primaries by Sparse Inversion via one-norm minimization
% This package is an application for simultaneous multiple removal and
% estimation of the source signature. It is based on the Estimation of Primaries
% by Sparse Inversion approach but reformulated via block- coordinate descent
% and convexification to L1 minimization. We believe this is therefore more
% robust to artifacts, and should converge quicker and more reliably that the
% original formulation.
%
% _This package is available only in the software release for members of SINBAD consortium._ 
%
% Tim Tai-Yi Lin (tlin@eos.ubc.ca),  April 30th 2013

%% How is the package organized?
% Robust EPSI is divided into two parts in this release.
%%
% The *program* itself (implemented in pure MATLAB) is contained in
% |{$SLIM_ROOT}/tools/algorithms/REPSI/|. However, its main driver file (|EPSI_SLIM_main.m|)
% is designed for experimentation, and requires a large myriad of options to be
% set. Even though only a few of the options are required, I nevertheless
% consider it too complicated to directly invoke as a utility.
%%
% Therefore, *interaction* with the Robust EPSI program should be done through
% the helper scripts in |{$SLIM_ROOT}/applications/WavefieldSeparation/RobustEPSI/|. Each
% subdirectory contained therein (except for |data| and |doc|) contain the
% helper scripts to necessary to process a specific dataset. See the file
% |README_examples.txt| for documentation on all the example datasets contained
% in this software release.


%% Getting up and running
% A basic walkthrough can be found <example.html here>.

%% Do I need the MATLAB Parallel Computation Toolbox?
% _No._ Robust EPSI is explicitly designed to work without the MATLAB Parallel
% Computation Toolbox. If your machine has enough memory to process the prestack
% data, it is perfectly fine to simply execute without parallel mode (setting
% |parallel = 0|). The total memory requirement is currently about 25x the
% stored single-precision filesize of SU data, or 80x if using the
% Curvelet-Wavelet-based transform domain.
%%
% Parallel execution mode for RobustEPSI mostly takes advantage of the
% shared-memory mode, and is meant to overcome the single-machine memory limit.
% If you lack the license necessary for distributed execution of parallel code,
% I recommend turning off parallel mode entirely to avoid the execution overhead
% associated with using distributed arrays. Because the Robust EPSI code is
% painstakingly vectorized, in most cases invoking parallel execution mode
% actually results in a slowdown for single-machine execution.

%% More information
% See the included |README.txt| and |README_example.txt| for a more detailed description of the program. 

%% References 
% <https://slim.gatech.edu/content/robust-estimation-primaries-sparse-inversion-one-norm-minimization [1]>. Tim T.Y. Lin and Felix J. Herrmann, [2013] Robust estimation of primaries by sparse inversion via one-norm minimization: Geophysics, *78*, no. 3, R133-R150
