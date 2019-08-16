%% set environment
curdir = pwd;
projdir = [pwd,'/../..'];
datadir = [projdir,'/data/'];
resultdir = [projdir,'/results/sigsbee/'];

%% i/o related
data_file = 'preproc_iwave_FDData.mat';
model_file = 'sigsbee_nosalt_model.mat';
result_file = 'RTMTotalFD_iwave_TrueQ.mat';

data_file = [datadir, data_file];
model_file = [datadir, model_file];
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

% S
S = eye(nsrc);
S = repmat(S,[1,1,nf]);
S_f = freq.^2.*exp(-(freq/freq_peak).^2).*exp(1i*2*pi* freq*(model.t0));
for i=1:nf
    S(:,:,i) = S(:,:,i)*S_f(i);
end
S = vec(distributed(S));

scale_op = 1/1000;
load(data_file,'P_up_SURF_t')
P_up_SURF_t = P_up_SURF_t/(-7.2);  % correct for scaling by
                                   % iwave modelling, compared to
                                   % SLIM Freq. domain modelling,
                                   % factor of 7.2 by LS matching
P_up_SURF = opFT*P_up_SURF_t(:);
P_up_SURF = reshape(P_up_SURF,[nf,nrec,nsrc]);
P_up_SURF = shiftdim(P_up_SURF,1);
P_up_SURF = vec(distributed(P_up_SURF));
Pr = scale_op*P_up_SURF;
mul_scale = 1/25;
Ps = -mul_scale*P_up_SURF;

clear P_up_SURF P_up_SURF_t

S = sourcegrid(nrec, nsrc, src_polarity, S);
Ps = sourcegrid(nrec, nsrc, src_polarity, Ps);
model.xsrc = model.xrec;
opBorn = oppDF_old(m0(:), S-Ps, model);
dm_RTM = opBorn'*Pr(:);
dm_RTM = reshape(dm_RTM,model.n);
dm_RTM = dm_RTM(model.np(1)+1:end,:);
save(result_file,'dm_RTM')
