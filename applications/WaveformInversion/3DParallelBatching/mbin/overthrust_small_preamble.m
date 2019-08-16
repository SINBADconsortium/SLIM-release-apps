%
% Author: Curt Da Silva, 2016
%
model_file = 'overthrust-small.rsf';

% Number of points to actually use
% 3.75 km x 3.75 km x 1km model
nx_s = 150; ny_s = 150; nz_s = 40;

entire_model = true;

% Free parameters
% Number of right hand sides over which to solve Helmholtz in
% parallel using multi-threading
numsrc_mt = 1;

% Number of sources (for x, y dimensions individually)
% Total number of sources is this number^2
numsrc = 10;

% Linear solve tolerance for FWI
linearsolve_tol = 1e-6;