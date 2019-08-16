function MADAGASCAR
slimcomp = getenv('SLIM_COMP');
if isempty(slimcomp)
	fprintf('FATAL ERROR: SLIM_COMP environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading MADAGASCAR Toolbox from\n\t%s\n',slimcomp);
if isdir([slimcomp '/external/lib'])
    assert(addSRpath('SLIM_COMP','/external/lib',false));
else
    warning('MATLAB-RSF is not present. Applications that require MATLAB-RSF will not work')
end
