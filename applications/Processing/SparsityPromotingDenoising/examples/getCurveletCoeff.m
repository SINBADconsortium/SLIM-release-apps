function getCurveletCoeff
% This script computes the curvelet coefficients of the input seismic data.
% Input data (this example): a frequency slice from the (2012) Gulf of Mexico data set

close all
mydir=pwd;

% Run the parameters script
setDenoiseParams;

% Load input data
if AcqGrid
   D = fdata2AcqGrid(fname_input, 1);
else
   D = rsf_read_all(fname_input);
end

% Data dimension
dim = size(D);

% Normalize the data 
D = D/nmfactor;

% Setup the measurement operator (the 2D curvelet operator)
% opCurvelet: (M, N, nbscales, nbangles, finest, ttype, is_real)
% Note: change parameters for opCurvelet here
C = opCurvelet(dim(1), dim(2), max(1,ceil(log2(min(dim(1),dim(2))) - 3)), nbangles, finest, ttype, is_real);

% Solve the one-norm recovery problem
% spgl1: (A, b, tau, sigma, x, options)
% xL1: estimated (synthesis) curvelet coefficients
xL1 = spgl1(C', D(:), opt.spgl1_tau, opt.spgl1_sigma, opt.spgl1_x, options); 

% Save curvelet coefficients
rsf_write_all(fname_curvcoeffs, {'out=stdout'}, xL1, [1 1])

cd(mydir);
