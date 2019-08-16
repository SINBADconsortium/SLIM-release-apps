function CurveLab
slimcomp = getenv('SLIM_COMP');
if isempty(slimcomp)
	fprintf('FATAL ERROR: SLIM_COMP environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	return;
end
fprintf('Loading CurveLab Toolbox from\n\t%s\n',slimcomp);
if isdir([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM'])
    assert(addSRpath('SLIM_COMP','/tools/transforms/CurveLab-2.1.2-SLIM/fdct_usfft_matlab',false));
    assert(addSRpath('SLIM_COMP','/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_matlab',false));
    assert(addSRpath('SLIM_COMP','/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_cpp/mex',false));
    assert(addSRpath('SLIM_COMP','/tools/transforms/CurveLab-2.1.2-SLIM/fdct3d/mex',false));
    assert(addSRpath('SLIM_COMP','/tools/transforms/CurveLab-2.1.2-SLIM/mecv',false));
else
    warning('CurvLab is not present. Applications that require CurveLab will not work')
end
