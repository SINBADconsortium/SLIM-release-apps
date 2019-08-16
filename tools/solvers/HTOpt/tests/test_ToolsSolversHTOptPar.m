if exist('runtests','file') == 0
    error('xunit package missing, please add it to the matlab path');
end
assert(parpool_size() > 0, 'Parallel pool is not on');

cdir=pwd;
tdir=fileparts(which('test_ToolsSolversHTOpt'));

cd(tdir);
runtests testHTuckOptPar;
cd(cdir);
