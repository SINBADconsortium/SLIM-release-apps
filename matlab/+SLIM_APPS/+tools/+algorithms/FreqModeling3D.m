function FreqModeling3D
slimapps = getenv('SLIM_APPS');
if isempty(slimapps)
	fprintf('FATAL ERROR: SLIM_APPS environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading 3DFreqModeling Toolbox from\n\t%s\n',slimapps);
assert(addSRpath('SLIM_APPS','/tools/solvers/Krylov',true));
assert(addSRpath('SLIM_APPS','/tools/solvers/Multigrid',false));
assert(addSRpath('SLIM_APPS','/tools/algorithms/CommonFreqModeling',false));
assert(addSRpath('SLIM_APPS','/tools/algorithms/3DFreqModeling',false));
