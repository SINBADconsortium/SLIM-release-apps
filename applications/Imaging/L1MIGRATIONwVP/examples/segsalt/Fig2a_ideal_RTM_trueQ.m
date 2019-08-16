% initialization and environment setup
curdir = pwd;
projdir = [pwd,'/../..'];
datadir = [projdir,'/data'];
resultdir = [projdir,'/results/segsalt'];

data_file = [datadir,'/data.mat'];
model_file = [datadir,'/model.mat'];
result_file = [resultdir,'/linear_RTM.mat'];

nt = 2001;
dt = 0.004;
df = 1/(nt*dt);
freq = df:df:12;

load(model_file)
m = vel_true;
m0 = vel_bg;
m = 1./(m.^2);
m0 = 1./(m0.^2);
dm = m-m0;

model.o = [0,0];
model.d = [24.384,24.384];
model.n = [160,645];
model.nb = [20,20];
model.freq = freq;
model.f0 = 5;
model.t0 = -0.25;
model.zsrc = model.d(1); 
model.xsrc = (0:2:(model.n(2)-1))*model.d(2);
model.zrec = model.d(1); 
model.xrec = (0:2:(model.n(2)-1))*model.d(2);

topmute_start = 3;
topmute_end = 6;

nsrc = length(model.xsrc);
nrec = length(model.xrec);
nf = length(model.freq);

Q = eye(nsrc);
opBorn = oppDF_old(m0(:),Q,model);

load(data_file,'D_linear');
D = D_linear;
clear D_linear

D = vec(distributed(reshape(D,nrec*nsrc,nf)));
dm_RTM = opBorn'*D;
dm_RTM = reshape(dm_RTM,model.n);

save(result_file,'dm_RTM')
