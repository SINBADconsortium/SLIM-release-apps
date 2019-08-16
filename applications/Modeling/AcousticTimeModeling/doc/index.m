%% acoustic time-domain modeling, linearized modeling and imaging
%
% The modeling code is based on a 9-point mixed-grid discretization of
% the 2D Helmholtz operator [1]. It solves the system in parallel over
% frequencies using direct factorization (Matlab's |mldivide|). Source
% injection and receiver sampling is done via cubic interpolation. The
% Jacobian is derived by linearizing the discretized system and its forward
% and adjoint action is calculated via the adjoint-state method.
%
% Author: Dr. Xiang Li
%
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
% the software release under |applications/Modeling/AcousticTimeModeling|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% The scripts can be run in serial mode but parallel mode is advised
% for the modeling and imaging examples. Use |parpool| to start parallel pool with the
% appropriate configuration and a divisor of 12 workers.

%% Examples
% A few examples are included here
%
% * Some examples of the modeling capabilities are shown in <TimeModeling_example.html TimeModeling_example.m>.
% * Basic tests of the modeling code are presented in <demo2D_gradient_test.html demo2D_gradient_test.m>


%% References
%	[1]. Xiang Li - 2D/3D time-stepping for imaging and inversion, SINBAD Fall consortium talks. 2014
