%% set environment
rng('default')
curdir = pwd;
projdir = [pwd,'/../..'];
datadir = [projdir,'/data/'];
resultdir = [projdir,'/results/sigsbee/'];

%% i/o related
data_file = 'EPSI3_iwave_FDData.mat';
model_file = 'sigsbee_nosalt_model.mat';
snapshot_file = 'MigTotalFD_Mul_iwave_26ss31freq_bg_GaussianEncoding_snapshots.mat';
result_file = 'MigTotalFD_Mul_iwave_26ss31freq_bg_GaussianEncoding.mat';

data_file = [datadir, data_file];
model_file = [datadir, model_file];
snapshot_file = [resultdir, snapshot_file];
result_file = [resultdir, result_file];

% time-frequency
nt = 1024;
dt = 0.008;
df = 1/(nt*dt);
freq_upper = 38;
freq_peak = 15;
freq = df:df:freq_upper;
freq_full = 0:df:df*floor(nt/2);
freq_mask = get_subset_mask(freq_full,freq);

% read model
load(model_file)
m0 = vel_bg;
m0 = 1e6./(m0.^2);
mute_start = 90;
mute_end = 100;

model.o = [0,0];
model.d = [7.62,7.62];
model.n = [500,781];
model.nb = [30 30];
model.freq = freq;
model.f0 = 0;
model.t0 = 0;
model.zsrc = model.d(1)*2; 
model.xsrc = (0:3:(model.n(2)-1))*model.d(2);
model.zrec = model.d(1)*0; 
model.xrec = (0:3:(model.n(2)-1))*model.d(2);
nsrc = length(model.xsrc);
nrec = length(model.xrec);
nf = length(model.freq);

% process for dipole source
src_polarity = 'di';
if strcmp(src_polarity,'di') 
    model.np = [model.zsrc/model.d(1)+1,0];
    model.o = [0,0]-model.d.*model.np;
    model.n = [500,781]+model.np;
    model.zsrc = model.zsrc*[-1, 1];
    m0 = [repmat(m0(1,:),[model.np(1),1]);m0];
    mute_start = mute_start+model.np(1);
    mute_end = mute_end+model.np(1);
else
    model.np = [0,0];
end

% Fourier transform operator
opFT1D = opFFTsym_conv_mask(nt,freq_mask);
opFT = opKron(opDirac(nrec*nsrc),opFT1D);

scale_op = 1/1000;
load(data_file,'multiples','initial_data')
Pr = multiples;
Pr = opFT*Pr(:);
Pr = reshape(Pr,[nf,nrec,nsrc]);
Pr = shiftdim(Pr,1);
Pr = vec(distributed(Pr));
Ps = initial_data;
Ps = opFT*Ps(:);
Ps = reshape(Ps,[nf,nrec,nsrc]);
Ps = shiftdim(Ps,1);
Ps = vec(distributed(Ps));
D.Pr = scale_op*Pr;
mul_scale = 1/25;
D.Ps = -mul_scale*Ps;
clear multiples initial_data Pr Ps

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
dep = dz:dz:dz*nz;
opDepth = opKron(opDirac(nx),opDiag(sqrt(dep)));
% top-muting op
opModelMute = opKron(opDirac(nx),opLinMute(nz,mute_start,mute_end));

iter = 50;
opts = spgSetParms('optTol',1e-3, ...
                   'bpTol' , 1e-3 , ...
                   'decTol',5e-2,...
                   'iterations',iter,...
                   'project', @NormL1_project, ...
                   'primal_norm', @NormL1_primal, ...
                   'dual_norm', @NormL1_dual,...
                   'quitPareto', 1,...
                   'linear',1,...
                   'lineSrchIt',2);
tau = 0;
sigma = 0;
x0 = zeros(size(C,2),1);

inv_para.nf_batchsize = 31;
inv_para.ns_batchsize = 26;
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

[x,info] = Mig_Mul(m0, D, model, inv_para);
dm = inv_para.opRight*x;
dm = reshape(dm, nz, nx);
dm = dm(model.np(1)+1:end,:);
save(result_file,'x','dm','info')
curvelet_denoise(result_file, 0.0015)