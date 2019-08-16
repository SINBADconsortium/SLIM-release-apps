function MADAGASCAR
slimcomp = getenv('SLIM_COMP');
if isempty(slimcomp)
	fprintf('FATAL ERROR: SLIM_COMP environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading NFFT Toolbox from\n\t%s\n',slimcomp);
if isdir([slimcomp '/external/share/nfft/matlab/nnfft'])
    assert(addSRpath('SLIM_COMP','/external/lib',false));
    assert(addSRpath('SLIM_COMP','/external/share/nfft/matlab/nnfft',false));
else
    warning('NFFT is not present. Applications that require NNFFT will not work')
end
