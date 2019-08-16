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
flagrandseq = 0;
flagsim = 1;
nsupshot = 2;

%% Working code


[bulkt n d o] = rsf_read_all(fname_bulkt);
bulkh        = rsf_read_all(fname_bulkh);
buoy         = rsf_read_all(fname_buoy);
dbulk        = rsf_read_all(fname_dbulk);
dbuoy        = rsf_read_all(fname_dbuoy);


m1    = [bulkt(:);buoy(:)];
m2    = [bulkh(:);buoy(:)];

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
%ldata      = ReadSuFast(fname_data,1001,'','b');
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
% Calling the wrapper function
D1     = Fd2DT(m1,model,options);
[D2,J] = Fd2DT(m2,model,options);

dm     = [dbulk(:);dbuoy(:)];

ddata  = J*dm;

n_data = [length(model.t), length(ddata(:))/length(model.t)];

figure;imagesc([1:n_data(2)],model.t,reshape(ddata,n_data));
title('Linearized Data','fontsize',16);
ylabel('Time [s]','fontsize',16);
colormap gray
colorbar

figure;imagesc([1:n_data(2)],model.t,reshape(D1-D2,n_data));
title('Non-Linearized Data','fontsize',16);
ylabel('Time [s]','fontsize',16);
colormap gray
colorbar

Dlin_rand  = reshape(ddata,n_data);
Dnlin_rand = reshape(D1-D2,n_data);

save('../result/D_lin_rand.mat','Dlin_rand');
save('../result/D_nlin_rand.mat','Dnlin_rand');




