%% test spot operator for iwave

%% Setting
% Toolbox path setting and data & model file setting
curdir = pwd;
%mbinpath = [curdir '/../../../../../tools/utilities/iWAVE'];
%addpath(mbinpath);
addpath(curdir);
% Important file names and their paths
fname_buoy  = [curdir '/../../demo_model/cambuoy.rsf'];
fname_bulkt = [curdir '/../../demo_model/cambulkt.rsf'];
fname_bulkh = [curdir '/../../demo_model/cambulkh.rsf'];
fname_dbulk = [curdir '/../../demo_model/camdbulk.rsf'];
fname_dbuoy = [curdir '/../../demo_model/camdbuoy.rsf'];

% General setting for using iWave++

options.fwdpara = 'fwdparams.txt';
options.linpara = 'linparams.txt';
options.delete   = 1;
options.bulkonly = 0;
options.nlab     = 4;

                        
%% Working code


[bulkt n d o] = rsf_read_all(fname_bulkt);
bulkh        = rsf_read_all(fname_bulkh);
buoy         = rsf_read_all(fname_buoy);
dbulk        = rsf_read_all(fname_dbulk);
dbuoy        = rsf_read_all(fname_dbuoy);


m1    = [bulkt(:);buoy(:)];
m2    = [bulkh(:);buoy(:)];

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
[D2,J] = Fd2DT(m2,model,options);

dm     = [dbulk(:);dbuoy(:)];

ddata  = J*dm;

figure;imagesc(model.xrec,model.t,reshape(ddata,length(model.t),length(model.xrec)));
title('Linearized Data','fontsize',16);
xlabel('Receiver [m]', 'fontsize', 16);
ylabel('Time [s]','fontsize',16);
colormap gray
colorbar

figure;imagesc(model.xrec,model.t,reshape(D1-D2,length(model.t),length(model.xrec)));
title('Non-Linearized Data','fontsize',16);
xlabel('Receiver [m]', 'fontsize', 16);
ylabel('Time [s]','fontsize',16);
colormap gray
colorbar
