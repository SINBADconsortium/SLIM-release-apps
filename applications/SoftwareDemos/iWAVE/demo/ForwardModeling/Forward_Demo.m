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
fname_bulkt = [curdir '/../../demo_model/cambulkt.rsf'];
%fname_data = '/users/slic/zfang/Project/IWAVE++/Orig/Rand_data.su';

% General setting for using iWave++

options.fwdpara  = 'fwdparams.txt';
options.delete   = 1;
options.nlab     = 4;



%% Working code

[bulkt n d o] = rsf_read_all(fname_bulkt);
buoy         = rsf_read_all(fname_buoy);



m1    = [bulkt(:);buoy(:)];

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
%ldata      = ReadSuFast(fname_data,1001,'','b');




% Calling the wrapper function
D1     = Fd2DT(m1,model,options);

figure;imagesc(model.xrec,model.t,D1);colorbar
title('Demo Shot Record','fontsize',16);
xlabel('Receiver [m]', 'fontsize', 16);
ylabel('Time [s]','fontsize',16);
colorbar 
colormap gray
