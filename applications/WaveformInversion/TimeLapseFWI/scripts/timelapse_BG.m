
clc;
addpath(genpath('./mbin/')); % functions needed to run the code

%setpath;

%% Load input data and model parameters
load ../data/ObsData
%% setup locations for results
expdir   = ['./results/' mfilename];
mkdir(expdir)
%% define initial velocity model
m = 1e6./(vbase(:)).^2;
S = opKron(opSmooth(n(2),50),opSmooth(n(1),50));
m0 = S*m(:);

vel_initial =  reshape(1e3./sqrt(m0),n);

%% Load/Set Inversion parameters
inv_params.beta = 1;%default is 1
inv_params.nproc = 16;%default is 8,16, depending on the number of processors in a node
inv_params.nsim = 10;% default is 10
inv_params.nexp = 10;% default is 10
inv_params.maxiter = 30; % maximum iterations for l1 solver, default is 30
             

%% Perform Inversion

[results] = FourD_FWI(vel_initial,Dbase,Dmon,model_base,model_monitor,inv_params,expdir);


%% Possible torque configuration
% j = torque.batchParallel('timelapse_BG',72,2,8,0);


