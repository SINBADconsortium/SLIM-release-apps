% Recover subsampled seismic lines by partitioning across offset slices
% and utilizing the correlations between the support of curvelet transform
% coefficients in the time-midpoint (TM) domain.
% 
% Copyright 2013 Hassan Mansour (hassanm@cs.ubc.ca)

close all;

% set working directories
label = 'wL1OffsetTM';
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
C   = opCurvelet(dim(1),dim(2), max(1,ceil(log2(min(dim(1),dim(2))) - 3)), 16, 1, 'ME', 1);

% subsample the time domain seismic line
for t = 1:dim(1)
    bH(t,:,:) = reshape(RM*D(t,:)',dim(2), dim(3));
end

display('Finished computing the measurements');

% build the mask volume
maskH_vol = permute(repmat(maskH, [1,1,dim(1)]),[3 1 2]);

%% interpolate the seismic line

% positive offsets
options.st = dim(2);
options.fin = dim(3);
options.partorder = [3 1 2];
options.transform = C;
options.maxiter = 500;
options.omega = 0.3;  % set omega = 1 to recover using standard L1 minimization.	

DHest = wL1min(bH, maskH_vol, options);


% negative offsets (or use symmetry when sources and receivers are aligned)
options.st = dim(2);
options.fin = 1;
options.partorder = [3 1 2];
options.transform = C;
options.maxiter = 500;
options.omega = 0.3;  % set omega = 1 to recover using standard L1 minimization.

[DHest(options.fin:options.st,:,:)] = wL1min(bH, maskH_vol, options);


% permute back to time, midpoint, offset
DHest = real(permute(DHest, [2 3 1]));

% convert back to the source receiver domain
for t = 1:dim(1) 
    Dest(t,:,:)=reshape(SR'*vec(real(DHest(t,:))),1,dimtime(2),dimtime(3));
end

%% save data in outputdir
D = Dtime;
save([outputdir, '/wL1_offset_TM.mat'], 'DH', 'DHest', 'maskH', 'D', 'Dest', 'mask', 'C', 'SR');

