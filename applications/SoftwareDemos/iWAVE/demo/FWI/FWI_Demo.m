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
fname_bulks = [curdir '/../../demo_model/cambulkh.rsf'];
fname_bulkt = [curdir '/../../demo_model/cambulkt.rsf'];
fname_data  = [curdir '/../../demo_model/cam_data.su'];
%fname_data = '/users/slic/zfang/Project/IWAVE++/Orig/Rand_data.su';


% General setting for using iWave++
options.fwdpara  = 'fwdparams.txt';
options.adjpara  = 'adjparams.txt';
options.delete   = 1;
options.bulkonly = 1;
options.nlab     = 16;




%% Working code

[bulks n d o] = rsf_read_all(fname_bulks);
buoy          = rsf_read_all(fname_buoy);
bulkt         = rsf_read_all(fname_bulkt);

m1    = [bulks(:);buoy(:)];
mt    = [bulkt(:);buoy(:)];

% Model Info

model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = 0:50:1000; 
model.zrec = 980;
model.xrec = 0:10:1000;
model.f0   =5;
model.t0   = .5 ;
model.n    = n;
model.o    = o;
model.d    = d;

%data       = ReadSuFast(fname_data,length(model.t),'','b');
data        = Fd2DT(mt,model,options);
data        = data(:);





% Setting for FWI
ncut    = 20;
fh           = @(x) misfit_cut(x,model,data(:),options,ncut);
opt.fid      = fopen(['lbfgs_iter.log'],'a');
opt.M        = 10;
opt.itermax  = 30;
opt.write    = 1;

% Use lbfgs to solve
mf           = lbfgs_wWolfe_modified(fh,m1,opt);             
bulkf        = reshape(mf(1:prod(model.n)),model.n);


figure;imagesc(model.n(2),model.n(1),bulkf);
title('FWI result','fontsize',16);
xlabel('X [m]','fontsize',16);
ylabel('Z [m]','fontsize',16);
colorbar

DFWI = bulkf;
%save('../result/D_FWI.mat','DFWI');


