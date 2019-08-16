% This script runs the RTM example. For parallel runs, use 3, 6 or 12
% workers. The script may take up to 2 hours to run.
%

setpaths;
curdir = pwd;
expdir = [resultsdir '/rtm'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end

% read models
[v,o,d,n]  = odnread([datadir '/marm_vp.odn']);

v = reshape(v,n);
v = [1500*ones(60,501); v(1:301,500 + [1:501])];
n = size(v);

[z,x] = odn2grid(o,d,n);

% model grid
model.o = o;
model.d = d;
model.n = n;

% absorbing boundary
model.nb = [50 50];

% source/receiver grid
model.zsrc = 25;
model.xsrc = 0:100:5000;
model.zrec = 15;
model.xrec = 0:50:5000;

% frequencies
model.freq = [2.5:.5:20];

% wavelet
model.f0 = 10;
model.t0 = 0;

% source matrix
Q = speye(length(model.xsrc));

% smoothing operators
S = opKron(opSmooth(n(2),200),opSmooth(n(1),200));

% smooth model
m  = 1e6./v(:).^2;
m0 = S*m;
dm = m - m0;

% make data
D       = F(m,Q,model);
[D0,J0] = F(m0,Q,model);
dD      = D - D0;
[od,dd,nd] = grid2odn(model.xrec,model.xsrc,model.freq);

% frequency weighting
Wf = oppKron2Lo(opDiag(1./model.freq),opDirac(nd(1)*nd(2)));

% depth weighting
Wz = opKron(opDirac(n(2)),opDiag(sqrt(z)));

% images

dmt = Wz*J0'*(Wf*dD);

odnwrite([expdir '/dmt.odn'],dmt,o,d,n);



