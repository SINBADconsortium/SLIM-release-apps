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

SLIM_APPS.tools.utilities.SPOT_SLIM;
SLIM_APPS.tools.utilities.pSPOT
SLIM_COMP.MADAGASCAR

% 2DFreqModeling
SLIM_APPS.tools.algorithms.FreqModeling2D;
SLIM_APPS.tools.solvers.Krylov;
SLIM_APPS.tools.solvers.QuasiNewton;
SLIM_APPS.tools.Miscellaneous;

% done
fprintf('Done loading SLIM Toolboxes\n');
%path
curdir = pwd;
addpath([curdir '/mbin']);

