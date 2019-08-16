% initialization and environment setup
rng('default')
curdir = pwd;
projdir = [pwd,'/../..'];
datadir = [projdir,'/data'];
resultdir = [projdir,'/results/segsalt'];

data_file = [datadir,'/data.mat'];
model_file = [datadir,'/model.mat'];
snapshot_file = [resultdir,'/linear_trueQ_GaussianEnc2_snapshots.mat'];
result_file = [resultdir,'/linear_trueQ_GaussianEnc2.mat'];

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

src_polarity = 'mo';
Q = eye(nsrc);
Q = repmat(Q,[1,1,nf]);
Q = vec(distributed(Q));

scale_op = 1/1000;
load(data_file,'D_linear');
D = vec(distributed(reshape(D_linear,[nrec,nsrc,nf])));
D = scale_op*D;
clear D_linear

% inversion
nz = model.n(1);
nx = model.n(2);
dz = model.d(1);
dx = model.d(2);
% curvelet synthesis op
NBscales = max(1,ceil(log2(min(nz,nx)) - 3));
NBangles = 16;
Finest = 1;
Ttype = 'ME';
IS_real = 0;
C = opCurvelet(nz,nx,NBscales,NBangles,Finest,Ttype,IS_real);
C = C';
% real restriction op
opGetReal = opRealRestriction(nx*nz);
% depth preconditioner
useDepthPrec = 1;
if useDepthPrec
    dep = dz:dz:dz*nz;
    opDepth = opKron(opDirac(nx),opDiag(sqrt(dep)));
else
    opDepth = opDirac(nx*nz);
end
% top-muting op
opModelMute = opKron(opDirac(nx),opLinMute(nz,topmute_start,topmute_end));

iter = 60;
opts = spgSetParms('iterations' , iter , ...
                   'verbosity' , 2 , ...
                   'bpTol' , 1e-3 , ...
                   'optTol' , 1e-3, ...
                   'quitPareto', 1,...
                   'decTol', 5e-2);
tau = 0;
sigma = 0;
x0 = zeros(size(C,2),1);

inv_para.nf_batchsize = 16;
inv_para.ns_batchsize = 16;
inv_para.f_renewal = 1;
inv_para.s_renewal = 1;
inv_para.encoding_type = 'gaussian';
inv_para.opRight = opModelMute*opDepth*opGetReal*C;
inv_para.opLeft = scale_op;
inv_para.tau = tau;
inv_para.sigma = sigma;
inv_para.x0 = x0;
inv_para.opts = opts;
inv_para.iter = iter;
inv_para.hub = 0.005;
inv_para.estS = 0;
inv_para.srmask = ones(nrec,nsrc);
inv_para.snapshot_file = snapshot_file;
inv_para.src_polarity = src_polarity;

[x,info] = Mig_Pri(m0, Q, D, model, inv_para);
dm = inv_para.opRight*x;
dm = reshape(dm, nz, nx);
save(result_file,'x','dm','info')

% curvelet denoising
op.nx = nx;
op.nz = nz;
op.dx = dx;
op.dz = dz;
op.nz_pad = 0;
op.mutestart = topmute_start;
op.muteend = topmute_end;
op.nbscales = NBscales;
op.nbangles = NBangles;
op.finest = Finest;
op.ttype = Ttype;
op.is_real = IS_real;
op.useDepthPrec = useDepthPrec;
curvelet_denoise_gen(result_file, 0.001, op)
