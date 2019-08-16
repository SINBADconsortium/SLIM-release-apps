% Recover interpolated baseline and monitor data (common-receiver gathers) and their 
% difference between  from a set of time-lapse data with randomly missing shots 
% between the gathers using the joint recovery method.


% Copyright 2014 Felix Oghenekohwo (foghenekohwo@eos.ubc.ca)
clc;
clear all;
close all;

% set working directories
label = 'MissShots_OneReceiverGather';
datadir = '../data/';
%toolsdir = '../../../../tools/algorithms/AdaptiveSparseRecovery';
outputdir = ['../results/' label];
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

% add paths to the working directories
%addpath(toolsdir);
addpath(datadir);

% load data and restrict to just one common-receiver gather
load([datadir, 'time-lapse-data.mat']);
D1 = squeeze(Baseline(:,75,:));
D2 = squeeze(Monitor(:,75,:));
clear Baseline Monitor

% build mask or choose a sampling mask applied in the source domain

dim = size(D1);

nrow = dim(1);
ncol = dim(2);

rho_1 = 0.5;
rho_2 = 0.5;

idx1 = randperm(ncol);
idx2 = randperm(ncol);
while idx2 == idx1
   idx2 = randperm(ncol);
   idx1 = randperm(ncol);
end


% define the mask for the independently selected randomly missing shots
I1 = opRestriction(ncol,idx1(1:ceil(rho_1*ncol)));
I2 = opRestriction(ncol,idx2(1:ceil(rho_2*ncol)));	
	
% define the independent sampling operators
R1 = opKron(I1,opDirac(nrow));
R2 = opKron(I2,opDirac(nrow));

% measure the observed data and combine both
y1 = R1*D1(:);
y2 = R2*D2(:);
y = [y1;y2];
Ry1 = R1'*y1;
Ry2 = R2'*y2;
% Build Operators
% curvelet 
C   = opCurvelet(nrow,ncol, max(1,ceil(log2(min(nrow,ncol)) - 3)), 16, 1, 'ME', 0);


% set up the JRM recovery operator and the synthesis operator
A  = [opStack(R1*C',R2*C')  opBlockDiag(R1*C',R2*C')];
S = [opStack(C',C')  opBlockDiag(C',C')];

% set options for spgl1 algorithm
options.iterations = 10;

% run the recovery algorithm
tic;xest_curv = spgl1(A, y, 0, 0, [], options);toc;

xr = real(S*xest_curv);

D1_rec = xr(1:nrow*ncol);
D2_rec = xr(nrow*ncol+1:end);

%% save data in outputdir
save([outputdir, '/SymmetricSamples.mat'], 'D1', 'D2', 'Ry1', 'Ry2', 'D1_rec', 'D2_rec','rho_1','rho_2');

