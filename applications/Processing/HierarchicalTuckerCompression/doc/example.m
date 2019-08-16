% Large-scale seismic data compression with on-the-fly shots/receivers generation from Hierarchical Tucker 
%
% This script displays the results after the expriment carried
% out in |applications/Processing/HierarchicalTuckerCompression/examples/BG_3D.m
% Users can change parameters to improve the results when doing compression
% or interpolaion in HT format
%
% Author: Yiming Zhang (yzhang@eoas.ubc.ca)
% Date: March 2018

baseDir    = pwd;
baseDir    = [baseDir(1:end-3)];
resultsDir = [baseDir 'results/'];
dataDir    = [baseDir 'data/'];

% Load previously computed results
load([resultsDir 'results.mat']);

% Load the data D with size (nrec x X nrec y X nsrc x X nsrc y)
load([dataDir '/BG_3Hz.mat']);

% No. of shots/receivers along x and y direction
nsrcx = size(D,3);
nrecx = size(D,1);
nsrcy = size(D,4);
nrecy = size(D,2);

% No. of shots/receivers in total
nrecs  =   nrecx * nrecy;
nsrcs  =   nsrcx * nsrcy;


%% Compression in fully sampled scenarios

% Permute the data into noncanical organization (receiver x, source x, receiver y, source y)
D     = permute(D,[1 3 2 4]);
% Compute the SNR of compressed data 
D1    = reshape(dimTree1.full(x1), [nrecx, nsrcx, nrecy, nsrcy]);
diff1 = D1 - D;
snr1  = SNR(D1,D);
disp(['SNR of entire compressed volume ' num2str(snr1) 'dB']);

% View the shot, shot1 is formed through on-the-fly shot extraction from HT
figure;
imagesc(reshape(real(shot1),[nrecx, nrecy]));colormap seiscol; caxis([-60 60]); colorbar
xlabel('receiver y'); ylabel('receiver x'); title('Compressed data')

figure;
imagesc(real(squeeze(D(:,6,:,2)))); colormap seiscol; caxis([-60 60]); colorbar
xlabel('receiver y'); ylabel('receiver x'); title('True data')

figure;
imagesc(real(squeeze(diff1(:,6,:,2)))); colormap seiscol; caxis([-60 60]); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('Difference')


%% Interpolation in missing entries scenarios 

% Compute the SNR of interpolated data 
D2    = reshape(dimTree2.full(x2),[nrecx, nsrcx, nrecy, nsrcy]);
diff2 = D2 - D;
snr2  = SNR(D2,D);
disp(['SNR of entire interpolated volume ' num2str(snr2) 'dB']);

% View the shot, shot2 is formed through on-the-fly shot extraction from HT
figure;
imagesc(reshape(real(shot2),[nrecx, nrecy]));colormap seiscol; caxis([-60 60]); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('Interpolated data')

figure;
imagesc(real(squeeze(D(:,6,:,2)))); colormap seiscol; caxis([-60 60]); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('True data')

figure;
imagesc(real(squeeze(diff2(:,6,:,2)))); colormap seiscol; caxis([-60 60]); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('Difference')


%% Compute D*w if you provide the probing vector 

% view the results of D*v1 from random vector v1
figure;
imagesc(reshape(real(d1), nrecx, nrecy)); colormap seiscol;caxis([-2 2]*1e3); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('From True D')
% view the results based on the idea of shots/receivers extraction from compressed HT after above compression, given the same v1
figure;
imagesc(reshape(real(d3), nrecx, nrecy)); colormap seiscol;caxis([-2 2]*1e3); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('From compressed HT')
% Diff
figure;
imagesc(reshape(real(d3 - d1), nrecx, nrecy)); colormap seiscol;caxis([-2 2]*1e3); colorbar 
xlabel('receiver y'); ylabel('receiver x'); title('Difference')

%% Compute D^H*w if you provide the probing vector 

% view the results of D^H*v2 from random vector v2
figure;
imagesc(reshape(real(d2), nsrcx, nsrcy)); colormap seiscol;caxis([-6 6]*1e3); colorbar 
xlabel('source y'); ylabel('source x'); title('From True D')
% view the results based on the idea of shots/receivers extraction from compressed HT after above compression, given the same v2
figure;
imagesc(reshape(real(d2), nsrcx, nsrcy)); colormap seiscol;caxis([-6 6]*1e3); colorbar 
xlabel('source y'); ylabel('source x'); title('From compressed HT')
% Diff
figure;
imagesc(reshape(real(d4 - d2), nsrcx, nsrcy)); colormap seiscol;caxis([-6 6]*1e3); colorbar 
xlabel('source y'); ylabel('source x'); title('Difference')


