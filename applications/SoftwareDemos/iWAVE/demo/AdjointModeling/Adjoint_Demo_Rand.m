%% test spot operator for iwave

%% Setting
% Toolbox path setting and data & model file setting
curdir = pwd;
%mbinpath = [curdir '/../../../../tools/utilities/iWAVE'];
%addpath(mbinpath);
addpath(curdir);
% Important file names and their paths
fname_buoy  = [curdir '/../../demo_model/cambuoy.rsf'];
fname_bulkh = [curdir '/../../demo_model/cambulkh.rsf'];
fname_bulkt = [curdir '/../../demo_model/cambulkt.rsf'];
fname_data  = [curdir '/../../demo_model/camd1_data.su'];
fname_dbulk = [curdir '/../../demo_model/camdbulk.rsf'];
fname_dbuoy = [curdir '/../../demo_model/camdbuoy.rsf'];

% General setting for using iWave++
options.adjpara  = 'adjparams.txt';
options.linpara  = 'linparams.txt';
options.delete   = 1;
options.bulkonly = 1;
options.nlab     = 4;
options.label    = 'Test';

flagrandseq = 0;
flagsim = 1;
nsupshot = 2;


%% Working code

[bulkt n d o] = rsf_read_all(fname_bulkt);
[bulkh n d o] = rsf_read_all(fname_bulkh);
buoy          = rsf_read_all(fname_buoy);
dbulk        = rsf_read_all(fname_dbulk);
dbuoy        = rsf_read_all(fname_dbuoy);


m1    = [bulkt(:);buoy(:)];
m2    = [bulkh(:);buoy(:)];

% Model Info

model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = [500 800];
model.zrec = 40;
model.xrec = 0:10:1000;
model.f0   = 30;
model.t0   = .2 ;
model.n    = n;
model.o    = o;
model.d    = d;

dt         = model.t(2)-model.t(1);
nt         = length(model.t);

if flagrandseq > 0
    model.randseq = 1;
end

if flagsim> 0
    for i = 1:nsupshot
        model.srcdata(:,:,i)   = vec(getrick(model.f0,model.t0,dt,nt)) * randn(1,length(model.xsrc));
    end
    model.simultsrc = 1;
end

J      = oppDFd2DT(m2,model,options);
dm     = [dbulk(:)];
%ldata         = ReadSuFast(fname_data,length(model.t),'','b');
ldata       = J * dm(:);


% Calling the wrapper function
%J      = oppDFd2DT(m,model,options);
gbulk  = J'*ldata(:);
gbulk  = reshape(gbulk,model.n);

figure;imagesc(model.n(2),model.n(1),gbulk);
title('Gradient','fontsize',16);
xlabel('X [m]','fontsize',16);
ylabel('Z [m]','fontsize',16);
colorbar
colormap gray

Dadj_rand = gbulk;
save('../result/D_adj_rand.mat','Dadj_rand');


