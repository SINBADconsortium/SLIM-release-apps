%% Large-scale seismic data compression with on-the-fly shots/receivers generation from compressed Hierarchical Tucker parameter
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software provides an algorithm for fast shots/receivers extraction
% once representing seismic frequency slices in the latent Hierarchical Tucker
% format. It can work for both fully sampled and missing entries scenarios.
% The combination of massive compression and fast on reconstruction of 3D shot
% or receiver gathers leads to a substantial reduction in memory costs in the
% subsequent shot-based processings.
%
% Author: Yiming Zhang (yzhang@eoas.ubc.ca)
%
% Date: March, 2018

%% Downloading & Dependencies
% The synthetic examples code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |applications/Processing/HierarchicalTuckerCompression|.
%
% The code has been tested with _Matlab R2015b_ and supports, but does not require the Parallel Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part of the 
% <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.

%% Running & Parallelism
% All of the examples and results are produced by the scripts found in this software release under |/applications/Processing/HierarchicalCompression/example|.
% Start matlab from |/applications/Processing/HierarchicalCompression| to add the appropriate paths.
%
% To run the scripts, follow the instructions in the README.md file enclosed with the code 

%% Examples and results
% <example.html Example is shown here.> You can reproduce the results by
% running the script
% under |/applications/Processing/HierarchicalCompression/examples|.

%% References
% <https://slim.gatech.edu/Publications/Public/Conferences/SEG/2017/zhang2017SEGmsd/zhang2017SEGmsd.pdf [1]> Y. Zhang, C. Da Silva, R. Kumar and F. J. Herrmann, 2017. Massive 3D seismic data compression and inversion with hierarchical Tucker
%
% <https://slim.gatech.edu/Publications/Private/Submitted/2018/silva2018alr/silva2018alr.pdf [2]> C. Da Silva, Y. Zhang, R. Kumar and F. J. Herrmann, 2018. Applications of low-rank compressed seismic data to full waveform inversion and extended image volumes

%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.