% Recover subsampled seismic lines by partitioning across frequency slices
% and utilizing the correlations between the support of curvelet transform
% coefficients in the midpoint-offset (MH) domain.
% 
% Copyright 2013 Hassan Mansour (hassanm@cs.ubc.ca)

close all;

% set working directories
label = 'wL1FreqMH';
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

% build mask or choose a sampling mask applied in the source-receiver
% domain
dim = size(D);
maskrecv = (randn(dim(2),1) >= 0);
mask = repmat(maskrecv, 1,dim(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% convert to midpoint-offset domain %%%%%%%%%%%%%%%%%
opSR=opSR2MH(dim(2));
info=opSR([],0);
SR=oppFunction(info{1},info{2},opSR);

maskH = logical(reshape(SR*int8(mask(:)),dim(2),2*dim(2)));

% Do conversion of full data to midpoint offset
for t = 1:dim(1) 
    DH(t,:,:)=reshape(opSR(real(D(t,:)),1),1,dim(2),2*dim(2));
end

% save the time domain data and dimension
Dtime = D; 
dimtime = size(Dtime); 

D = DH;
dim = size(D);

% set up the sampling operator in the MH domain
RM=opMask(dim(2)*dim(3),find(maskH));

% Build Operators
% curvelet 
C   = opCurvelet(dim(2),dim(3), max(1,ceil(log2(min(dim(2),dim(3))) - 3)), 16, 1, 'ME', 1);
% 1D Fouerier
F = opDFT(dim(1));

% build the frequency seismic line
Df = F*D;

display('Finished converting time axis to frequency axis');


% want to recover only the zero and positive frequencies since data is real
Dftrunk = Df(1:floor(dim(1)/2)+1,:,:);

dimf = size(Dftrunk);
% subsample the frequency seismic line
for w = 1:dimf(1)
    bf(w,:,:) = reshape(RM*Dftrunk(w,:)',dimf(2), dimf(3));
end

display('Finished calculating the measurements');

% build the mask volume
maskH_vol = permute(repmat(maskH, [1,1,dim(1)]),[3 1 2]);


%% interpolate the seismic line
options.st = 1;
options.fin = dimf(1);
options.transform = C;
options.maxiter = 500;
options.omega = 0.3;  % set omega = 1 to recover using standard L1 minimization.

[Dfest] = wL1min(bf, maskH_vol, options);

%% recover the negative frequencies by taking the conjugate symmetric flip
% of the positive frequencies
Dfestflip = flipdim(conj(Dfest(2:end, :,:)),1);
Dfestfull = [Dfest(1:end-1,:,:);Dfestflip];

display('converting frequency seismic line to time seismic line');

% convert back to the time domain
DHest = real(F'*conj(Dfestfull));

% convert back to the source receiver domain
for t = 1:dim(1) 
    Dest(t,:,:)=reshape(SR'*vec(real(DHest(t,:))),1,dimtime(2),dimtime(3));
end

%% save data in outputdir
D = Dtime;
save([outputdir, '/wL1_freq_MH.mat'], 'DH', 'DHest', 'maskH', 'D', 'Dest', 'mask', 'C', 'SR');
