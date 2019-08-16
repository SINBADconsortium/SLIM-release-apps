function pSPOT
slimapps = getenv('SLIM_APPS');
if isempty(slimapps)
	fprintf('FATAL ERROR: SLIM_APPS environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading pSPOT Toolbox from\n\t%s\n',slimapps);
assert(addSRpath('SLIM_APPS','/tools/utilities/pSPOT',false));
try
   isdir(fullfile(spot.path,'tests','xunit'));
catch ME
   error('Can''t find SPOT or its xunit toolbox.');
end
