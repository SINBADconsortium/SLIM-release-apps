% Plotting results from using JRM
% Copyright 2014 Felix Oghenekohwo (foghenekohwo@eos.ubc.ca)

clear all;
close all;

% set working directories
label = 'MissShots_OneReceiverGather';
setpath;
outputdir = [resultsdir '/' label];

if ~exist(outputdir,'dir')
    error('JRM_OneRecvGather.m  must be run before plotting figures.');
end

load([outputdir, '/SymmetricSamples.mat']);

dim = size(D1);
D1_rec = reshape(D1_rec,dim(1),dim(2));
D2_rec = reshape(D2_rec,dim(1),dim(2));

Ry1 = reshape(Ry1,dim(1),dim(2));
Ry2 = reshape(Ry2,dim(1),dim(2));


cax = 0.5;
cmap = 'gray';

% Original data 
figure; 
subplot(1,3,1);imagesc(D1); colormap(cmap); caxis([-cax cax]); 
title('Baseline'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,2);imagesc(D2); colormap(cmap); caxis([-cax cax]); 
title('Monitor'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,3);imagesc(D1-D2); colormap(cmap); caxis(0.1*[-cax cax]); 
title('Difference'); xlabel('Source number'); ylabel('Time samples');


% Measured data 
figure; 
subplot(1,2,1);imagesc(Ry1); colormap(cmap); caxis([-cax cax]); 
title('Baseline'); xlabel('Source number'); ylabel('Time samples');
subplot(1,2,2);imagesc(Ry2); colormap(cmap); caxis([-cax cax]); 
title('Monitor'); xlabel('Source number'); ylabel('Time samples');
%subplot(1,3,3);imagesc(Ry1-Ry2); colormap(cmap); caxis(0.1*[-cax cax]); 
%title('Difference'); xlabel('Source number'); ylabel('Time samples');


% Recovered data 
figure; 
subplot(1,3,1);imagesc(D1_rec); colormap(cmap); caxis([-cax cax]); 
title('Baseline'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,2);imagesc(D2_rec); colormap(cmap); caxis([-cax cax]); 
title('Monitor'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,3);imagesc(D1_rec-D2_rec); colormap(cmap); caxis(0.1*[-cax cax]); 
title('Difference'); xlabel('Source number'); ylabel('Time samples');


% Difference
figure; 
subplot(1,3,1);imagesc(D1-D1_rec); colormap(cmap); caxis([-cax cax]); 
title('Baseline'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,2);imagesc(D2-D2_rec); colormap(cmap); caxis([-cax cax]); 
title('Monitor'); xlabel('Source number'); ylabel('Time samples');
subplot(1,3,3);imagesc((D1-D2)-(D1_rec-D2_rec)); colormap(cmap); caxis(0.1*[-cax cax]); 
title('Difference'); xlabel('Source number'); ylabel('Time samples');