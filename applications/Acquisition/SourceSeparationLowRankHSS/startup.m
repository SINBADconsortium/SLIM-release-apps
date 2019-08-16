%%% start SLIM_COMP

% MADAGASCAR
SLIM_COMP.MADAGASCAR


%%% start SLIM_APPS

% Miscellaneous
SLIM_APPS.tools.Miscellaneous

% Low-rank minimization
SLIM_APPS.tools.algorithms.LowRankMinimization

% Generalized SPGL1
SLIM_APPS.tools.solvers.GenSPGL1

% SPOT (slim updates)
SLIM_APPS.tools.utilities.SPOT_SLIM

% pSPOT
SLIM_APPS.tools.utilities.pSPOT

% SegyMAT
SLIM_APPS.tools.utilities.segyMAT

%%% done
fprintf('Done loading SLIM Toolboxes\n');


% Private path
curdir = pwd;
funcsdir = [curdir, '/misc_funcs'];
addpath(funcsdir)

