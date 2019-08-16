%% test spot operator for iwave
clear all;

%% Setting
% Toolbox path setting and data & model file setting
curdir = pwd;
%mbinpath = [curdir '/../../../../../tools/utilities/iWAVE'];
%addpath(mbinpath);
addpath(curdir);
% Important file names and their paths
fname_buoy  = [curdir '/../../demo_model/cambuoy.rsf'];
fname_bulkh = [curdir '/../../demo_model/cambulkh.rsf'];
fname_data  = [curdir '/../../demo_model/camd1_data.su'];

% General setting for using iWave++
options.adjpara  = 'adjparams.txt';
options.delete   = 1;
options.bulkonly = 1;
options.nlab     = 4;



%% Working code


[bulkh n d o] = rsf_read_all(fname_bulkh);
buoy          = rsf_read_all(fname_buoy);


m    = [bulkh(:);buoy(:)];

% Model Info

model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = 500;
model.zrec = 40;
model.xrec = 0:10:1000;
model.f0   = 30;
model.t0   = .2 ;
model.n    = n;
model.o    = o;
model.d    = d;

ldata         = ReadSuFast(fname_data,length(model.t),'','b');



% Calling the wrapper function
J      = oppDFd2DT(m,model,options);
gbulk  = J'*ldata(:);
gbulk  = reshape(gbulk,model.n);

figure;imagesc(model.n(2),model.n(1),gbulk);
title('Gradient','fontsize',16);
xlabel('X [m]','fontsize',16);
ylabel('Z [m]','fontsize',16);
colorbar
colormap gray

