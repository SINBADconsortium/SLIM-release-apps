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


%%% start SLIM_APPS
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
% WRI
addpath([slimapps '/tools/algorithms/WRI']);
% Miscellaneous
addpath([slimapps '/tools/operators/misc']);
addpath([slimapps '/tools/functions/misc']);
addpath([slimapps '/tools/functions/misc/ndim_dist']);

% my extras
%addpath([slimapps '/applications/WaveformInversion/2DWRI/mbin']);
addpath(genpath([slimapps '/applications/WaveformInversion/2DWRI-SourceEstimation/mbin']));
addpath([slimapps '/applications/WaveformInversion/ConstrainedFWI/mbin']);
addpath(genpath([slimapps '/tools/solvers/QuasiNewton/minConf_mod']))
%addpath([slimapps '/tools/algorithms/CommonFreqModeling']);
%addpath([slimapps '/tools/algorithms/2DFreqModeling']);
addpath(genpath([slimapps '/tools/solvers/LeastSquares']));
%addpath([slimapps '/tools/algorithms/3DFreqModeling']);
addpath([slimapps '/tools/solvers/Krylov']);
%%% done
fprintf('Done loading SLIM Toolboxes\n');
%path
