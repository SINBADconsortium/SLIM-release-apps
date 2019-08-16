fprintf('\n****************************');
fprintf('\n****** STARTING SERIAL TESTS\n\n');

fprintf('\n ***** SPOT tests start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.utilities.SPOT_SLIM
spottests;

fprintf('\n ***** SPGL1 test start in 5 seconds\n\n');
clear;
pause(5);
pause off;
SLIM_APPS.tools.solvers.SPGL1_SLIM
spgdemo;
pause on;

fprintf('\n ***** Quasi-Newton tests start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.QuasiNewton
test_ToolsSolversQuasiNewton;

fprintf('\n ***** Tests in tools/solvers/HTOpt/tests start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.HTopt
test_ToolsSolversHTOpt;

fprintf('\n ***** Tests in tools/solvers/LinearizedBregman start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.LinearizedBregman
test_ToolsSolversLinearizedBregman;

fprintf('\n ***** Tests in tools/algorithms/2DFreqModeling start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.Krylov
SLIM_APPS.tools.algorithms.FreqModeling2D
test_ToolsAlgorithms2DFreqModeling;

fprintf('\n ***** Tests in tools/algorithms/3DFreqModeling start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.Krylov
SLIM_APPS.tools.algorithms.FreqModeling3D
test_ToolsAlgorithms3DFreqModeling;

fprintf('\n ***** Tests in tools/algorithms/TimeModeling start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.algorithms.TimeModeling
test_TimeModeling_serial;

fprintf('\n ***** Tests in tools/operators/misc start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.Miscellaneous
test_ToolsOperatorsMisc;

SLIM_COMP.CurveLab
if numel(which('fdct_sizes_mex'))
    fprintf('\n ***** 2D Curvelet tests start in 5 seconds\n\n');
    clear;
    pause(5);
    A=opCurvelet(64,64);
    disp(A.utest());
else
    fprintf('\n WARRNING! ***** 2D Curvelet NOT INSTALLED - skipping test\n\n');
end
if numel(which('fdct3d_sizes_mex'))
fprintf('\n ***** 3D Curvelet tests start in 5 seconds\n\n');
    clear;
    pause(5);
    A=opCurvelet3d(64,64,64);
    disp(A.utest());
else
    fprintf('\n WARRNING! ***** 3D Curvelet NOT INSTALLED - skipping test\n\n');
end

SLIM_COMP.MADAGASCAR
if numel(which('rsf_write_all'))
    fprintf('\n ***** MATLAB-RSF tests start in 5 seconds\n\n');
    clear;
    pause(5);
    X=magic(100);
    rsf_write_all('test.rsf',{},X);
    Y=rsf_read_all('test.rsf');
    !sfrm test.rsf
    rsf_write_all('test.rsf',{'out=stdout'},X);
    Z=rsf_read_all('test.rsf');
    !sfrm test.rsf
    if norm(X-Y)==0 && norm(X-Z)==0
	    fprintf('PASSED!\n');
    else
	    fprintf('FAILED RSF READ/WRITE test!\n');
    end
else
    fprintf('\n WARRNING! ***** MATLAB-RSF NOT INSTALLED - skipping test\n\n');
end

fprintf('\n****** All serial tests done');
fprintf('\n****************************\n');
exit;
