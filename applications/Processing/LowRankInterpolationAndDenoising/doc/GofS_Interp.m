%% Seismic data regularization, missing-trace interpolation and denoising
%
% This script display the output results of regularization, missing-trace interpolation
% and denoising.
% Please note that the regularization, interpolation and denoised results are produced by 
% running the example scripts without using the optimal parameters. 
% Users can play with parameter to improve the results.
%
% Author: Rajiv Kumar (rakumar@eos.ubc.ca)
% Date: April, 2014

%% Data dimensions
curdir     = pwd;
basedir    = curdir(1:end-4);
datadir    = [basedir '/data/'];
resultsdir = [basedir '/results/'];

% Data parameter
ntime = 1024;
nrec  = 355;
nsrc  = 355;
d     = [0.004 12.5 12.5];
t     = (0:1:(ntime-1))*d(1);
xrec  = (0:1:(nrec-1))*d(2);
xsrc  = (0:1:(nsrc-1))*d(3);
idx = 177;
%% Regularization
% Perform regularization and missing-trace interpolation by running the example file in the directory 
% |/applications/Processing/LowRankInterpolationAndDenoising/examples/GulfofSuez/GofS_regularize_and_interp.m|
options.datafile = [datadir 'Interpolation/SuezShots125.355shots.su'];
% Read input data
D              = ReadSuFast(options.datafile);
D              = reshape(D,ntime,nrec,nsrc);

% Read output results after regularization 
load([resultsdir 'GulfofSuez/GulfofSuezreginterp.mat']);
Dt = reshape(Dt,ntime,nrec+1,nsrc+1);
% Display regularization results
% View a shot gather
figure;
imagesc(xrec,t,squeeze(D(:,:,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Ground truth');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,squeeze(Dt(:,1:nrec,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Regularized shot gather');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,(squeeze(D(:,:,idx))-squeeze(Dt(:,1:nrec,idx))));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Residual');colormap('gray');caxis([-1 1]*2e2)
%% Missing Trace Interpolation
% Perform interpolation by running the example file in the directory 
% |/applications/Processing/LowRankInterpolationAndDenoising/examples/GulfofSuez/GofS_Interp.m|

% Read output results after interpolation
load([resultsdir 'GulfofSuez/GulfofSuezinterp.mat']);

% Display interpolated results
% View a shot gather
figure;
imagesc(xrec,t,squeeze(D(:,:,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Ground truth');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,squeeze(output(:,:,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Interpolated shot gather');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,(squeeze(D(:,:,idx))-squeeze(output(:,:,idx))));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Residual');colormap('gray');caxis([-1 1]*2e2)

%% Missing Trace Interpolation and denoising
% Perform interpolation by running the example file in the directory 
% |/applications/Processing/LowRankInterpolationAndDenoising/examples/GulfofSuez/GofS_Interp_and_denoise.m|

% Read output results after Interpolation
load([resultsdir 'GulfofSuez/GulfofSuezinterpanddenoise.mat']);

% Display interpolated and denoised results
% View a shot gather
figure;
imagesc(xrec,t,squeeze(D(:,:,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Ground truth');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,squeeze(output(:,:,idx)));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Interpolated and denoised shot gather');colormap('gray');caxis([-1 1]*2e2)
figure;
imagesc(xrec,t,(squeeze(D(:,:,idx))-squeeze(output(:,:,idx))));
xlabel('Receiver(m)');ylabel('time(s)');ylim([0 2]);
title('Residual');colormap('gray');caxis([-1 1]*2e2)
