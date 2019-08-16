% read SLIM_ROOT environment
slimcomp = getenv('SLIM_COMP');
if length(slimcomp)<1
	fprintf('FATAL ERROR: SLIM_COMP environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
slimapps = getenv('SLIM_APPS');
if length(slimapps)<1
fprintf('FATAL ERROR: SLIM_APPS environment does not exist\n');
fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
return;
end

% start
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimapps);

% SPOT (slim updates)
addpath([slimapps '/tools/utilities/SPOT-SLIM']);
addpath([slimapps '/tools/utilities/SPOT-SLIM/tests']);
try
addpath(fullfile(spot.path,'tests','xunit'));
catch ME
error('Can''t find xunit toolbox.');
end

% pSPOT
addpath([slimapps '/tools/utilities/pSPOT']);

% CurveLab
if isdir([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM'])
	addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_usfft_matlab']);
	addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_matlab']);
	addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_cpp/mex']);
	addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct3d/mex']);
	addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/mecv']);
else
    %warning('CurveLab directory not found');
	error('CurveLab directory not found');
end
% SPGL1-SLIM (parallel)
addpath([slimapps '/tools/solvers/SPGL1-SLIM']);

% Quasi-Newton ()
addpath([slimapps '/tools/solvers/QuasiNewton']);

% Krylov Solvers
addpath([slimapps '/tools/solvers/Krylov']);

% CommonFreqModeling
addpath([slimapps '/tools/algorithms/CommonFreqModeling']);

% 2DFreqModeling
addpath([slimapps '/tools/algorithms/2DFreqModeling']);

% 3DFreqModeling
addpath([slimapps '/tools/algorithms/3DFreqModeling']);

% REPSI
addpath([slimapps '/tools/algorithms/REPSI']);

% Adaptive Sparse Recovery
addpath([slimapps '/tools/algorithms/AdaptiveSparseRecovery']);

% Low Rank Minimization
addpath([slimapps '/tools/algorithms/LowRankMinimization']);

% Miscellaneous
addpath([slimapps '/tools/operators/misc']);
addpath([slimapps '/tools/functions/misc']);

% segyMAT
addpath([slimapps '/tools/utilities/SegyMAT']);

% MADAGASCAR
if isdir([slimcomp '/external/lib'])
	addpath([slimcomp '/external/lib']);
else
    %warning('MADAGASCAR directory not found');
	error('MADAGASCAR directory not found');
end
% done
fprintf('Done loading SLIM Toolboxes\n');

% private path
curdir = pwd;
acqfuncsdir = [curdir, '/acq_funcs'];
addpath(acqfuncsdir)

