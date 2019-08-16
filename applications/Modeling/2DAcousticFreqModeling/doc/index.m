%% 2D constant-density acoustic frequency-domain modeling, linearized modeling and imaging
%
% The modeling code is based on a 9-point mixed-grid discretization of
% the 2D Helmholtz operator [1]. It solves the system in parallel over
% frequencies using direct factorization (Matlab's |mldivide|). Source
% injection and receiver sampling is done via cubic interpolation. The
% Jacobian is derived by linearizing the discretized system and its forward
% and adjoint action is calculated via the adjoint-state method.
%
%

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
% _algorithms/2DFreqModeling_
%
% * |Helm2D| - Construct Helmholtz matrix
% * |F| - modeling operator 
% * |DF| - Jacobian
% * |G| - modeling operator using analytic solution for constant and linear velocity profiles
% * |legendreQ| - evaluate Legendre Q function (used for |G|).
%
% _operators/misc_
%
% * |opLInterp1D| - 1D cubic Lagrange interpolation
% * |opLInter2D|  - 2D linear Lagrange interpolation
% * |opExtension| - Pads input with zeros or constant values
% * |opSmooth| - 1D smoothing by convolution with triangular kernel
% * |opSpline1D| - 1D cubic spline evaluation
% * |opSpline2D| - 2D cubic spline evaluation
% * |opDFTR|     - FFT for real input, outputs only positive frequencies.
%
% _functions/misc_
% 
% * |grid2odn|, |odn2grid| - convert grid vectors to [origin, increment, size] triplet and vice
% versa
% * |vec|, |invvec| - vectorize multidimensional array and reshape vector
% into multidimensional array.

%% Examples
% A few examples are included here
%
% * Some examples of the modeling capabilities are shown in <modeling.html modeling.m>.
% * Basic tests of the modeling code are presented in <testing.html testing.m>
% * Some imaging examples are found in <imaging.html imaging.m>


%% References
% <http://dx.doi.org/10.1190/1.1443979 [1]> C-H Jo,* C. Shin,* and J.H. Suh, 1996. An optimal 9-point, finite-difference, frequency-space, 2-D scalar wave extrapolator
% Geophysics 61(2), 529-537.
