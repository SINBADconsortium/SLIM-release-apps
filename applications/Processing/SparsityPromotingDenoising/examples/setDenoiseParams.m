% This script sets the parameters for denoising. 

% Set paths
setpath;
cd(resultsdir);

% Input data file
fname_input = [datadir '/FreqInd30_DG.rsf'];  

% Determine acquisition geometry
% AcqGrid: 1 - static geometry (SG) => fixed receivers
%          0 - dynamic geometry (DG) => moving receivers
AcqGrid = 0;

% Parameters for data
% Normalize the input data if it has very low amplitude
nmfactor = 1e-6;

% Set parameters for the 2D curvelet operator
% opCurvelet: (M, N, nbscales, nbangles, finest, ttype, is_real)
% Note - M, N and nbscales depend on the input data (see the script: getCurveletCoeff.m)
%      - type of curvelet transform can be either 'WRAP' (wrapping curvelets) or 'ME' (mirror-extended curvelets). See the README file.
nbangles = 16;
finest = 1;
ttype = 'WRAP';
is_real = 0;

% Set a label for the experiment
if strcmp(ttype,'WRAP')
   if AcqGrid
      label = 'FreqInd30_SG_WRAPcurv';
   else
      label = 'FreqInd30_DG_WRAPcurv';
   end
elseif strcmp(ttype,'ME')
   if AcqGrid
      label = 'FreqInd30_SG_MEcurv';
   else
      label = 'FreqInd30_DG_MEcurv';
   end
end

% Set options for the SPGL1 solver
% see function: spgl1(A, b, tau, sigma, x, options)
opt.spgl1_tau = 0;
opt.spgl1_sigma = 0;
opt.spgl1_x = [];
options.fid = fopen([label '.log'], 'w');
options.verbosity = 1;
options.iterations = 200;
options.optTol = 1e-04;
options.ignorePErr = 1;

% Output file names
fname_curvcoeffs = [label '_curvcoeffs.rsf'];
fname_denoised = [label '_denoised.rsf'];

% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

