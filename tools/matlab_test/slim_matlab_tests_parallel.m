fprintf('\n******************************');
fprintf('\n****** Starting parallel tests\n\n');

try
	pool=parpool_open(2);
catch mepool
	disp(mepool.message)
	fprintf('FATAL ERROR: unable to open parallel MATLAB pool\n');
	fprintf('             Tests did not run. Check the info above\n');
	fprintf('\n****** All paralle tests skipped');
	fprintf('\n********************************\n');
	exit;
end

fprintf('\n ***** pSPOT tests start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.utilities.pSPOT
pSPOT.utils.disableUselessWarnings(1)
pspottests;
pSPOT.utils.disableUselessWarnings(0)

fprintf('\n ***** Tests in tools/solvers/HTOpt/tests start in 5 seconds\n\n');
clear;
pause(5);
SLIM_APPS.tools.solvers.HTopt
test_ToolsSolversHTOptPar;

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
test_TimeStepping_parallel;

fprintf('\n');
parpool_close();
fprintf('\n****** All parallel tests done');
fprintf('\n******************************\n');
exit;
