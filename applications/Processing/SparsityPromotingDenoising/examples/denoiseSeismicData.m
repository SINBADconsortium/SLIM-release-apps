function denoiseSeismicData
% This script runs the denoising algorithm.

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
C = opCurvelet(dim(1), dim(2), max(1,ceil(log2(min(dim(1),dim(2))) - 3)), nbangles, finest, ttype, is_real);

% Load the curvelet coefficients
xL1 = rsf_read_all(fname_curvcoeffs);

% Denoising via debiasing 
% Note: check the function 'selectCoeffSupp_denoise' to see how different values of perc affect the denoised result. Change the values of perc until the denoised result looks reasonable. 
perc = 0.9922;  % This means 99.22%
[AT, Dden] = selectCoeffSupp_denoise(C, xL1, perc, D, nmfactor);

% If AcgGrid = 1, rearrange the denoised result on to dynamic grid 
% Note: check the function 'fdata2AcqGrid' to see how the variables change according to the input
if AcqGrid
   mode = 2;
   ncorig = 4001;
   startIND = 801;
   nr_subset = 801;
   nc_subset = 3201;
   Dden = fdata2AcqGrid(Dden, mode, ncorig, startIND, nr_subset, nc_subset);
end

% Reload input data for plotting
clear D
D = rsf_read_all(fname_input);

% View results
% Plotting parameters
cax = 1e-6;
cmap = seiscol;
xlab = 'Shot (#)';
ylab = 'Channel (#)';

% Input data
figure; imagesc(real(D), [-1 1]*cax); colormap(cmap); title('Input data')
xlabel(xlab); ylabel(ylab);

% Denoised data
figure; imagesc(real(Dden), [-1 1]*cax); colormap(cmap); title('Denoised data')
xlabel(xlab); ylabel(ylab);

% Difference data
figure; imagesc(real(D - Dden), [-1 1]*cax); colormap(cmap); title('Difference')
xlabel(xlab); ylabel(ylab);


% Save denoised data
rsf_write_all(fname_denoised, {'out=stdout'}, Dden, [1 1], [0 0], {'Channel' 'Shot'}, {'#' '#'})

cd(mydir);
