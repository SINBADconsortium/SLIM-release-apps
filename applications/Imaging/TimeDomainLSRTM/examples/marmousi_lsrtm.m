%% Sparsity promoting LSRTM with linearized Bregman
%
% Example script for running least-squares RTM with sparsity promotion on 
% the Marmousi model. 
%
% Author: Philipp Witte
% Date: February 2016
%
%

%% Set paths to result,data,model etc.

rng('default')

% paths  and name of results and snapshots
snapshotX = 'UpdateTest_X.mat';
snapshotZ = 'UpdateTest_Z.mat';
result = 'Result.mat';
params = 'Parameters.mat';

% pull data and model
if ~exist('../data/modelMarmousi.mat','file')
	cd ../data/; ! scons; cd ../examples/;
end

% Marmousi slowness^2 model [s^2/km^2]
load ../data/modelMarmousi.mat

% Load 2D data cube [ns x nrec x nsrc]
load ../data/linDataMarmousi.mat

%% Modeling parameters

% model parameters
model.o = [0, 0, 0];  % Origin [m]
model.n = [320, 800, 1];  % no. of grid points [z,x,y]
model.d = [10, 10, 10];   % grid spacing [m]
model.freesurface = 0;

% data parameters
model.f0   = 0.030; % peak frequency [MHz]
model.T= 4000;  % Recording length [ms]
model.NyqT=0:4:model.T; % Shot record time axis [ms]

% source/receiver position [m]
model.xsrc = 25:25:(n(2)*d(2));
model.zsrc = 10*ones(size(model.xsrc));
model.ysrc =  0*ones(size(model.xsrc));
model.xrec = 5:10:n(2)*d(2)-5;
model.zrec = 210;	
model.yrec = 0;	

% Code specific options
model.save = 'RAM'; % save wavefields in memory
model.type='full'; % acquisition geometry 
model.space_order=4; % Order of Laplacian

% Domain decomposition
model.ddcompx = 1;
model.ddcompz = 1;
model.ddcompy = 1;

%% LSRTM parameters

% LSRTM options
options.iterations = 40;
options.lambda = [];
options.sigma = 0;
options.numShots = 8;   % No. of shots per Bregman iteration
options.snapshotX = snapshotX;
options.snapshotZ = snapshotZ;

% Preconditioners
precon.modelDepth = 1;
precon.modelTopmute = 1;
precon.dataScaling = 1;
precon.dataTopmute = 1;


%% Setup grid and modeling options

% Calculate optimal grid size and time step
[m,model,m0]=Setup_CFL(m,model,m0);

% Source wavelet
q=sp_RickerWavelet(model.f0,1/model.f0,model.dt,model.T);


%% Call algorithm and save results
[dm,info] = linbregLSRTM(dData,model,m0,q,precon,options);

save(result,'dm','info');
save(params,'model');




