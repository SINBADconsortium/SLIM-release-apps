
%% test spot operator for iwave
clear all;

%% Setting
% Toolbox path setting and data & model file setting
curdir = pwd;
%mbinpath = [curdir '/../../../../../tools/utilities/iWAVE'];
%addpath(mbinpath);
addpath(curdir);
addpath('/home/slim/zfang/Project/SofterRelease/SLIM-release-apps/tools/utilities/iWAVE');
% Important file names and their paths
fname_buoy  = [curdir '/../../demo_model/layerbuoy.rsf'];
fname_bulks = [curdir '/../../demo_model/layerbulks.rsf'];
fname_data  = [curdir '/../../demo_model/layer_data.su'];
%fname_data = '/users/slic/zfang/Project/IWAVE++/Orig/Rand_data.su';

% General setting for using iWave++
options.fwdpara  = 'fwdparams.txt';
options.adjpara  = 'adjparams.txt';
options.delete   = 0;
options.bulkonly = 1;
options.nlab     = 16;

               
                        
%% Working code


[bulks n d o] = rsf_read_all(fname_bulks);
buoy         = rsf_read_all(fname_buoy);


m1    = [bulks(:);buoy(:)];

% Model Info

model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = 0:50:1000; 
model.zrec = 40;
model.xrec = 0:10:1000;
model.f0   = 30;
model.t0   = .2 ;
model.n    = n;
model.o    = o;
model.d    = d;

data       = ReadSuFast(fname_data,length(model.t),'','b');

% Setting the parameters, these will override the values of these variables



% Calling the wrapper function
[D1,J]     = Fd2DT(m1,model,options);

gbulk      = J'*(data(:)-D1(:));

gbulk      = reshape(gbulk,model.n);

figure;imagesc(model.n(2),model.n(1),gbulk);
title('RTM Imaging','fontsize',16);
xlabel('X [m]','fontsize',16);
ylabel('Z [m]','fontsize',16);
colorbar
colormap gray

