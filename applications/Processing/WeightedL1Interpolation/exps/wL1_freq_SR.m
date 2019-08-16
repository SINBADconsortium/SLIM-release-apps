% Recover subsampled seismic lines by partitioning across frequency slices
% and utilizing the correlations between the support of curvelet transform
% coefficients in the source-receiver (SR) domain.
% 
% Copyright 2013 Hassan Mansour (hassanm@cs.ubc.ca)

close all;

% set working directories
label = 'wL1FreqSR';
datadir = '../data/';
toolsdir = '../../../../tools/algorithms/AdaptiveSparseRecovery';
outputdir = ['../results/' label];
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

% add paths to the working directories
addpath(toolsdir);
addpath(datadir);


% load data and restrict to the first 500 time samples
load([datadir, 'GulfOfSuez.mat']);
D = D(1:500,:,:);

% build mask
dim = size(D);
maskrecv = (randn(dim(2),1) >= 0);
mask = repmat(maskrecv, 1,dim(3));

RM=opMask(dim(2)*dim(3),find(mask));

% Build Operators
% curvelet 
C   = opCurvelet(dim(2),dim(3), max(1,ceil(log2(min(dim(2),dim(3))) - 3)), 16, 1, 'ME', 1);
% 1D Fouerier
F = opDFT(dim(1));

% build the frequency seismic line
Df = F*D;

display('Finished converting time axis to frequency axis');

% want to recover only the zero and positive frequencies
Dftrunk = Df(1:floor(dim(1)/2)+1,:,:);

dimf = size(Dftrunk);
% subsample the frequency seismic line
for w = 1:dimf(1)
    bf(w,:,:) = reshape(RM*Dftrunk(w,:)',dimf(2), dimf(3));
end

display('Finished calculating the measurements');

% build the mask volume
mask_vol = permute(repmat(mask, [1,1,dim(1)]),[3 1 2]);

%% interpolate the seismic line
options.st = 1;
options.fin = dimf(1);
options.transform = C;
options.maxiter = 500;
options.omega = 0.3;  % set omega = 1 to recover using standard L1 minimization.

[Dfest] = wL1min(bf, mask_vol, options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % x = F*f;
% % x2 = x(1:end/2+1);
% % x3 = flipud(conj(x2(2:end)));
% % xtilde = [x2(1:end-1);x3];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 
Dfestflip = flipdim(conj(Dfest(2:end, :,:)),1);
Dfestfull = [Dfest(1:end-1,:,:);Dfestflip];

display('converting frequency seismic line to time seismic line');

% convert back to the time domain
Dest = real(F'*conj(Dfestfull));

%% save data in outputdir
save([outputdir, '/wL1_freq_SR.mat'], 'D', 'Dest', 'mask', 'C');
