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


flagrandseq = 0;
flagsim = 1;


nsupshot = 2;

% General setting for using iWave++

options.fwdpara  = 'fwdparams.txt';
options.delete   = 1;
options.nlab     = 4;







%% Working code

[bulkt n d o] = rsf_read_all(fname_bulkt);
buoy         = rsf_read_all(fname_buoy);



m1    = [bulkt(:);buoy(:)];

% Model Info
model   = [];
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
dt = model.t(2) - model.t(1);
nt = length(model.t);

if flagrandseq > 0
    model.randseq = 1;
    D1     = Fd2DT(m1,model,options);
end

if flagsim> 0
    for i = 1:nsupshot
        model.srcdata(:,:,i)   = vec(getrick(model.f0,model.t0,dt,nt)) * randn(1,length(model.xsrc));
    end
    model.simultsrc = 1;
    D1     = Fd2DT(m1,model,options);
end




%figure;imagesc(model.xrec,model.t,D1);colorbar
figure;imagesc(D1);colorbar
title('Demo Shot Record','fontsize',16);
xlabel('Receiver [m]', 'fontsize', 16);
ylabel('Time [s]','fontsize',16);
colorbar 
colormap gray

Dfwd_rand = D1;
save('../result/D_fwd_rand.mat','Dfwd_rand');
