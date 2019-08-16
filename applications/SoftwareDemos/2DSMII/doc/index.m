%% A parallel matrix-free framework for frequency-domain seismic modelling, imaging and inversion in Matlab
%
% These webpages proved a demonstration of the modeling, imaging and
% inversion framework. For more details please refer to the accompanying
% paper:
%
% T. van Leeuwen, 2012. <https://slim.gatech.edu/node/27141 A parallel matrix-free framework for frequency-domain seismic modelling, imaging and inversion in Matlab>.
%
%% Download
% This is part of the Public release of the SLIM Software,
% for download instructions, see <https://slim.gatech.edu/releases>
%
%% Dependencies
% The modeling code uses the following packages, found in the |tools| part
% of the software release.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.

%% Running & Parallelism
% All the examples can be reproduced by running the scripts found in
% the software release under |applications/SoftwareDemos/2DSMII|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% The scripts can be run in serial mode but parallel mode is advised
% for the modeling and imaging examples. Use |parpool| to open parallel pool with the
% appropriate configuration.

%% Functions
% The framwork consists of 3 distinct packages which can be found in the |tools| part of
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
% * |odnread|, |odnwrite| - read and write for gridded data.
%
% _functions/solvers/QuasiNewton_
%
% * |lbfgs| - LBFGS method with a weak Wolfe linesearch.
%
%

%% Examples
% The examples presented in the paper are found here.
%
% * Some tests and examples of the modeling capabilities are shown in <modeling.html modeling.m>.
% * Basic tests of the modeling code are presented in <testing.html testing.m>
% * Some inversion examples are found in <examples.html examples.m>
