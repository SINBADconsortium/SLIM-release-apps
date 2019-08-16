if exist('runtests','file') == 0
    error('xunit package missing, please add it to the matlab path');
end

cdir=pwd;
tdir=fileparts(which('test_ToolsSolversHTOpt'));

cd(tdir);
runtests testHTuck;
runtests testHTuckOpt;
cd(cdir);
