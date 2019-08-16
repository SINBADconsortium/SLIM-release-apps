%% 2D/3D acoustic anisotropic time-domain modeling and linearized modeling
%
% The modeling operator is based on a star 1D stencil of order 2,4 or 6. 
% It solves the system in parallel over sources . 
% Source injection and receiver sampling is done via cubic interpolation and
% exponential damping over a 3x3 square around the source location. The
% Jacobian is derived by linearizing the discretized system and its forward.
%
% Author: Mathias Louboutin, Philipp Witte 
%
%         September 2015 
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia

%% Dependencies
% The modeling code uses the following packages, found in the |tools| part
% of the software release.
%
% * _SPOT_ - object oriented framework for matrix-free linear algebra.
% * _pSPOT_ - parallel extension of SPOT.

%% Running & Parallelism
% All the examples can be reproduced by running the scripts found in
% the software release under |applications/Modeling/2DAcousticFreqModeling|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% The scripts can be run in serial mode but parallel mode is advised
% for the modeling and imaging examples. Use |parpool| to start parallel pool with the
% appropriate configuration and a divisor of 12 workers.

%% Functions
% The modeling code consists of 3 distinct packages which can be found in the |tools| part of
% the software release. The main components are listed below
%
% _algorithms/TimeModeling_
%
% * |opF| - modeling operator 
% * |opJ| - Jacobian
% * |Gen_data| - Data generation function
% * |Born| - Born modelling and RTM function
% * |GS| - function that output [f,g] misfit and gradient for FWI
%
% _operators/misc_
%
% * |opLInterp1D| - 1D cubic Lagrange interpolation
% * |opExtension| - Pads input with zeros or constant values
% * |opSmooth| - 1D smoothing by convolution with triangular kernel
%
% _functions/misc_
% 
% * |grid2odn|, |odn2grid| - convert grid vectors to [origin, increment, size] triplet and vice
% versa

%% Examples
% A few examples are included here
%
% * Some examples of the 2D modeling capabilities are shown in <Modelling_TTI_2D.html Modelling_TTI_2D.m>.
% * Some examples of the 3D modeling capabilities are shown in <Modelling_3D.html Modelling_3D.m>


%% References
% <https://slim.gatech.edu/Publications/Private/TechReport/2015/witte2015TRoam/witte2015TRoam.html .[1]> Philipp Witte*, Mathias Louboutin and Felix J. Herrmann, Overview on anisotropic modeling and inversion, Technical report
