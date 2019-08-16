%%% read SLIM_COMP/APPS environments
slimcomp = getenv('SLIM_COMP');
if length(slimcomp)<1
        fprintf('FATAL ERROR: SLIM_COMP environment does not exist\n');
        fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
        return;
end
slimapps = getenv('SLIM_APPS');
if length(slimapps)<1
        fprintf('FATAL ERROR: SLIM_APPS environment does not exist\n');
        fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
        return;
end


%%% start SLIM_COMP
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimcomp);

% CurveLab
if isdir([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM'])
    addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_usfft_matlab']);
    addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_matlab']);
    addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_cpp/mex']);
    addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/fdct3d/mex']);
    addpath([slimcomp '/tools/transforms/CurveLab-2.1.2-SLIM/mecv']);
else
    warning('CurvLab is not present. Applications that require CurveLab will not work')
end
% MADAGASCAR
if isdir([slimcomp '/external/lib'])
    addpath([slimcomp '/external/lib']);
else
    warning('MATLAB-RSF is not present. Applications that require MATLAB-RSF will not work')
end


%%% start SLIM_APPS
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimapps);

% SPOT (slim updates)
addpath([slimapps '/tools/utilities/SPOT-SLIM']);
addpath([slimapps '/tools/utilities/SPOT-SLIM/tests']);
try
   addpath(fullfile(spot.path,'tests','xunit'));
catch ME
   error('Can''t find xunit toolbox.');
end
% pSPOT
addpath([slimapps '/tools/utilities/pSPOT']);
% SPGL1-SLIM (parallel)
addpath([slimapps '/tools/solvers/SPGL1-SLIM']);
% Generalized SPGL1 ()
addpath([slimapps '/tools/solvers/GenSPGL1']);
addpath([slimapps '/tools/solvers/GenSPGL1/Penalty']);
% Quasi-Newton ()
addpath([slimapps '/tools/solvers/QuasiNewton']);
% HTopt package
addpath([slimapps '/tools/solvers/HTOpt']);
addpath([slimapps '/tools/solvers/HTOpt/htuck']);
addpath([slimapps '/tools/solvers/HTOpt/htuckOpt']);
addpath([slimapps '/tools/solvers/HTOpt/tuckOpt']);
addpath([slimapps '/tools/solvers/HTOpt/interpolation']);
addpath([slimapps '/tools/solvers/HTOpt/utility']);
addpath([slimapps '/tools/solvers/HTOpt/objective']);
% Krylov + Multigrid Solvers
addpath(genpath([slimapps '/tools/solvers/Krylov']));
addpath(genpath([slimapps '/tools/solvers/Multigrid']));
% segyMAT
addpath([slimapps '/tools/utilities/SegyMAT']);
% CommonFreqModeling
addpath([slimapps '/tools/algorithms/CommonFreqModeling']);
% 2DFreqModeling
addpath([slimapps '/tools/algorithms/2DFreqModeling']);
% 3DFreqModeling
addpath(genpath([slimapps '/tools/algorithms/3DFreqModeling']));
addpath(genpath([slimapps '/applications/Modeling/3DAcousticFreqModeling/scripts']));
% WRI
addpath([slimapps '/tools/algorithms/WRI']);
% REPSI
addpath([slimapps '/tools/algorithms/REPSI']);
% Adaptive Sparse Recovery
addpath([slimapps '/tools/algorithms/AdaptiveSparseRecovery']);
% Low Rank Minimization
addpath([slimapps '/tools/algorithms/LowRankMinimization']);
% Miscellaneous
addpath([slimapps '/tools/operators/misc']);
addpath([slimapps '/tools/functions/misc']);
%  addpath([slimapps '/tools/functions/linear_algebra']);


%%% done
fprintf('Done loading SLIM Toolboxes\n');
%path 
