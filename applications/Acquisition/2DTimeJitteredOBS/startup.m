%%% start SLIM_COMP

% CurveLab
SLIM_COMP.CurveLab

% MADAGASCAR
SLIM_COMP.MADAGASCAR


%%% start SLIM_APPS

% Miscellaneous
SLIM_APPS.tools.Miscellaneous

% SPGL1-SLIM (parallel)
SLIM_APPS.tools.solvers.SPGL1_SLIM

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
acqfuncsdir = [curdir, '/acq_funcs'];
addpath(acqfuncsdir)

