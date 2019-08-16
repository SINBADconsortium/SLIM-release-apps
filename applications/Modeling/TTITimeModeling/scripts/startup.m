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


%%% start SLIM_COMP
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimcomp);

%%% start SLIM_APPS
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimapps);

% SPOT (slim updates)
SLIM_APPS.tools.utilities.SPOT_SLIM
% pSPOT
SLIM_APPS.tools.utilities.pSPOT
% TimeModeling
SLIM_APPS.tools.algorithms.TimeModeling
% Misc
SLIM_APPS.tools.Miscellaneous

%%% done
fprintf('Done loading SLIM Toolboxes\n');
%path
