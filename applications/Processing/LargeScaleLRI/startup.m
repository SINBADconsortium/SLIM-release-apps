%%% read SLIM_COMP/APPS environments
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
% SPGLR_PAR package
addpath([slimapps '/tools/solvers/SPGLR_PAR']);
% Experiment parameters
addpath([slimapps '/applications/Processing/LargeScaleLRI/spgLR_experiments']);
% Miscellaneous
addpath([slimapps '/tools/operators/misc']);
addpath([slimapps '/tools/functions/misc']);


%%% done
fprintf('Done loading SLIM Toolboxes\n');
%path
