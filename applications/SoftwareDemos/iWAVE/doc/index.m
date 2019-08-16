%% Matlab Interface for iWAVE++
%
% This software release includes demos to use iWAVE to simultate time domain acoustic data, born modeling data, adjoint modeling,
% RTM, least-square migration and FWI in Matlab
%
% [1] Symes, W.W., Sun, D. and Enriquez, M., 2011. From modelling to inversion: designing a well‐adapted simulator. Geophysical Prospecting, 59(5), pp.814-833.
%
% [2] Tschannen,V., Fang, Z., and Herrmann, F. J., 2014. Time domain least squares migration and dimensionality reduction, Technical report, UBC.
%
% Author: Zhilong Fang
% Date  : Feb, 2016


%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/SoftwareDemos/iWAVE|.
%
% The code has been tested with _Matlab R2014b_ and requires the iWAVE in the Madagascar version 9520 and installation of extra software from SLIM-release-comp repository.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _utilities/iWAVE_ - toolbox for using iWAVE
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/SoftwareDemos/iWAVE/|. 
% Start matlab from that directory.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% Most scripts can be run in serial mode and parallel mode. When run the parallel mode, you do not need to open matlab pool
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Demos
% The demos specific to this package can be found in the |demo| directory.
% 
% * |_ForwardModeling/Forward_Demo|       - Demo script to simulate the acoustic wave data 
% * |_LinearizedModeling/Linear_Demo|     - Demo script to simulate the linearized data
% * |_AdjointModeling/Adjoint_Demo|       - Demo script to perform the adjoint operator of the linearized operator
% * |_RTM/RTM_Demo|                       - Demo script to run a simple RTM
% * |_LSM/LSM_Demo|                       - Demo script to run a simple least square migration
% * |_FWI/FWI_Demo|                       - Demo script to run a simple FWI



%% Results 
% <./example.html The results are shown here>.


