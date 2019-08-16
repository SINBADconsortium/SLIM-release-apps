%% Efficient least-squares imaging with sparsity promotion and compressive sensing
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes an parallel framework in Matlab for L1 
% migration [1,2], which based on the ideas from compressive-sensing and
% stochastic optimization, where the least-squares imaging result are 
% computed from random subsets of the data via curvelet-domain 
% sparsity-promotion. There are two different subset sampling strategies 
% are considered in this package: randomized source encoding, and drawing 
% sequential shots firing at random source locations from marine data 
% with missing near and far offsets. In both cases, we obtain excellent 
% inversion results compared to conventional methods at reduced 
% computational costs. There is also a small example based on a small one-
% layer model which can allow users to test the algorithm in a short time 
%
% Author: Xiang Li
% Date  : April, 2013

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware/ SLIM sofware release> under
% |/applications/Imaging/L1MIGRATIONbasic/|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware/ SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _transforms/CurveLab-2.1.2-SLIM/_ - Curvelet transform toolbox


%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/Imaging/L1MIGRATIONbasic/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% Scripts to reproduce imaging from one-layer model,  as well as
% results from sevaral papers [1],[2] can be foundin the |scripts| directory. <examples.html
% The results are shown here>.
%

%% References
%
% <https://slim.gatech.edu/content/efficient-least-squares-imaging-sparsity-promotion-and-compressive-sensing [1]> Felix J. Herrmann and Xiang Li, “Efficient least-squares imaging with sparsity promotion and compressive sensing”, Geophysical Prospecting, vol. 60, p. 696-712, 2012. 
%
% <https://slim.gatech.edu/content/randomized-dimensionality-reduction-full-waveform-inversion [2]> Felix J. Herrmann and Xiang Li, “Randomized dimensionality reduction for full-waveform inversion”, in EAGE Technical Program Expanded Abstracts, 2010.
%
%% Acknowledgements
% The synthetic Compass model was provided by the BG-GROUP, see also the
% <BG_DISCLAIMER.txt disclaimer>.
