%% Matlab interface for iWAVE++
%
%% Forward modeling

curdir = pwd;

mtfile = '../demo_model/cambulkt.rsf'; % model file
dfile  = '../result/D_fwd.mat';       % data file
[m n d o] = ReadAllData(mtfile);
t = 0:.004:2;
xrec = 0:10:1000;
load(dfile);
[z x] = odn2grid(o,d,n);
figure;imagesc(x,z,m);colorbar;xlabel('x [m]');ylabel('z [m]');title('bulk modulus model')
figure;imagesc(xrec,t,Dfwd);colorbar;xlabel('receiver [m]');ylabel('time [m]');title('Data');colormap redblue;
colorbar;caxis([-1,1]*1e2);

%% Linearize modeling

fname_dbulk = [curdir '/../demo_model/camdbulk.rsf'];
fname_bulkt = [curdir '/../demo_model/cambulkt.rsf'];
fname_bulkh = [curdir '/../demo_model/cambulkh.rsf'];
dfile1      = ['../result/D_lin.mat'];
dfile2      = ['../result/D_nlin.mat'];

[mt n d o] = ReadAllData(fname_bulkt);
[mh n d o] = ReadAllData(fname_bulkh);
[dm n d o] = ReadAllData(fname_dbulk);

load(dfile1);
load(dfile2);
[z x] = odn2grid(o,d,n);
figure;imagesc(x,z,mt);colorbar;xlabel('x [m]');ylabel('z [m]');title('true bulk modulus model');
figure;imagesc(x,z,mh);colorbar;xlabel('x [m]');ylabel('z [m]');title('background bulk modulus model');
figure;imagesc(x,z,dm);colorbar;xlabel('x [m]');ylabel('z [m]');title('bulk modulus difference');
figure;imagesc(xrec,t,Dnlin);colorbar;xlabel('receiver [m]');ylabel('time [m]');title('nonlinear data difference');colormap redblue;caxis([-15 15])
figure;imagesc(xrec,t,Dlin);colorbar;xlabel('receiver [m]');ylabel('time [m]');title('linear data difference');colormap redblue;caxis([-15 15]);

%% Adjoint
dfile = ['../result/D_adj.mat'];
load(dfile);
figure;imagesc(x,z,Dadj);colorbar;xlabel('x [m]');ylabel('z [m]');title('adjoint operator result');colormap redblue;caxis([-300 300]);

%% RTM
fname_bulks = [curdir '/../demo_model/layerbulks.rsf'];
fname_bulkt = [curdir '/../demo_model/layerbulkt.rsf'];
dfile = ['../result/D_RTM.mat'];
[mt n d o] = ReadAllData(fname_bulkt);
[ms n d o] = ReadAllData(fname_bulks);
load(dfile);
[z x] = odn2grid(o,d,n);
figure;imagesc(x,z,mt);colorbar;xlabel('x [m]');ylabel('z [m]');title('true bulk modulus model');
figure;imagesc(x,z,ms);colorbar;xlabel('x [m]');ylabel('z [m]');title('background bulk modulus model');
figure;imagesc(x,z,mt-ms);colorbar;xlabel('x [m]');ylabel('z [m]');title('model difference');colormap redblue
figure;imagesc(x,z,DRTM);colorbar;xlabel('x [m]');ylabel('z [m]');title('RTM result');colormap redblue;caxis([-1,1]*1e4);


%% LSM
dfile = ['../result/D_LSM.mat'];
load(dfile);
figure;imagesc(x,z,DLSM);colorbar;xlabel('x [m]');ylabel('z [m]');title('Least square migration result');colormap redblue;caxis([-1,1]*1e-1);

%% FWI
dfile = ['../result/D_FWI.mat'];
fname_bulkt = [curdir '/../demo_model/cambulkt.rsf'];
fname_bulks = [curdir '/../demo_model/cambulkh.rsf'];
[mt n d o] = ReadAllData(fname_bulkt);
[ms n d o] = ReadAllData(fname_bulks);
load(dfile);
figure;imagesc(x,z,mt);colorbar;xlabel('x [m]');ylabel('z [m]');title('true bulk modulus model');a=caxis;
figure;imagesc(x,z,ms);colorbar;xlabel('x [m]');ylabel('z [m]');title('initial bulk modulus model');caxis(a);
figure;imagesc(x,z,DFWI);colorbar;xlabel('x [m]');ylabel('z [m]');title('FWI result');caxis(a);








