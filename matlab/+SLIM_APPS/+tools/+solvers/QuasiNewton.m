function QuasiNewton
slimapps = getenv('SLIM_APPS');
if isempty(slimapps)
	fprintf('FATAL ERROR: SLIM_APPS environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading Quasi-Newton Toolbox from\n\t%s\n',slimapps);
assert(addSRpath('SLIM_APPS','/tools/solvers/QuasiNewton',true));
