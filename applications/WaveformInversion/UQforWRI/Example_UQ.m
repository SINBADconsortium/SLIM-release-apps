%% Uncertainty quantification (UQ) for Wavefield Reconstruction Inversion (WRI)
%This is a synthetic example of UQ for WRI . We mimic the reflection case.
%Sources and receivers are both loacted at the top of the model.


% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: April, 2018.

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (zfang@eoas.ubc.ca)

%This script will generate data and invert using the same modeling code.
%Model is the Camembert.

expk = 1;
curdir = pwd;

% Create the Result directory
expdir               = ['./result'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end

% Set up the file path for data and models
vpfile               = [curdir '/data/Ini_Vel.mat'];
vtfile               = [curdir '/data/SimpleLayerkm.mat'];
datafile             = [curdir '/data/SimpleLayerData.mat'];
sigmafile            = [curdir '/data/SigmaD.mat'];

% Set up the penalty parameters
lambda    = ones(9,1)*1e6; % 2.2731e+04 2.0900e+04 1.9540e+04
lambda(4) = 2.2731e+03;
lambda(5) = 2.0900e+03;
lambda(6) = 1.9540e+03;

% Set up the prior covariance matrix
beta      = 0;
sigma_D   = ReadAllData(sigmafile);
sigma_D   = sigma_D(1) * ones(9,1);
sigma_p   = 1;
itermax   = 50;
PriorType = 'smooth';
nb        = 20;
Priormp         = 1;
flagPriorVel    = 1;
nlayercst       = 0;
velprior        = 1;
nsmp            = 1e3;

[velp n d o]    = ReadAllData(vpfile);
[veltr]         = ReadAllData(vtfile);
vel             = velp;
velp            = velp(:);
m               = 1./vel(:).^2;
SIGMA_mp        = ones(n) * 0.33;
SIGMA_mp        = SIGMA_mp(:);
a               = min(SIGMA_mp.^2) * 0.9;
c               = SIGMA_mp.^2 - a;
b               = 0.65;
dmax            = 2;
C               = SmoothCov2(n,a,b,c,dmax);
R               = chol(C);
Lpm             = opInverse(R');
Lpm2            = opInverse(R);
ISIGMAp         = Lpm2*Lpm;


% Set up the source function and frequency band
f0       = 6;
t0       = .2;
If       = {[4  :  6]};
k0       = 1;

Cfile   = mfilename('fullpath');
Cfile   = [Cfile '.m'];
copyfile(Cfile,expdir);
cd(expdir)

[Dobs nd dd od] = ReadAllData(datafile);
Do              = reshape(Dobs,nd(2),nd(4),nd(5));


% setup parameters
model.o  = o;
model.d  = d;
model.n  = n;
model.nb = nb*[1 1];
[model.zrec,model.xrec,model.zsrc,model.xsrc,model.freq] = odn2grid(od,dd,nd);
%W            = speye(prod(model.n));
model.f0        = f0;
model.t0        = t0;
model.sigma     = sigma_D;
model.sigma_p   = sigma_p;
model.PriorType = PriorType;
model.flagPriorVel = flagPriorVel;
model.Priormp   = Priormp;
model.beta      = beta;
model.mp        = velp(:);
model.SIGMA_mp  = SIGMA_mp;
model.velprior  = velprior;
model.SIGMAp    = ISIGMAp;
Q               = speye(prod(nd(3:4)));
opts.itermax    = itermax;
opts.write      = 1;
opts.M          = 10;

params.wri      = 1;
params.hessian  = 'sparse';
params.nthreads = 6;
params.dist_mode = 'freq';
params.lambda    = lambda;
params.nsmp      = 8;
params.alpha_CI  = .95;
params.tol       = 1e-6;
params.maxiter   = 1000;
params.distribute = 0;
params.extend     = 0;
params.tmp_name   = ['TMP_SAMP_' num2str(expk)];



fun_df      = @(x) -2*x(:).^(-3);
Hvs         = diag(fun_df(vel(:)));

veli = (veltr(:) + velp(:))/2;
m    = 1./(veli(:).^2);

% Compute the MAP estimate
for k = k0:length(If)
    opts.fid         = fopen(['iter' num2str(k) '.log'],'w');;
    modelk           = model;
    paramsk          = params;
    modelk.freq      = model.freq(If{k});
    modelk.sigma     = model.sigma(If{k});
    modelk.SIGMA_mp  = modelk.SIGMA_mp * length(If);
    paramsk.lambda   = params.lambda(If{k});
    Dok              = Do(:,:,If{k});
    Dok              = distributed(Dok);
    Dok              = Dok(:);
    fh               = @(m) misfit_uq_slowness2(m,Dok,Q,modelk,paramsk);
    m                = lbfgs_wWolfe(fh,m,opts);
    mkdir(['./xm' num2str(k)]);
    movefile('*.dat',['./xm' num2str(k)]);
    velf             = reshape(1./sqrt(m),model.n);
    WriteAllData(['m_' num2str(k) '.mat'], velf, n, d, o);
end

WriteAllData(['MAP.mat'], velf, n, d, o);

params.hessian  = 'hess-app';

% Draw samples
for k = k0:length(If)
    opts.fid         = fopen(['iter' num2str(k) '.log'],'w');
    modelk           = model;
    paramsk          = params;
    modelk.freq      = model.freq(If{k});
    modelk.sigma     = model.sigma(If{k});
    paramsk.lambda   = params.lambda(If{k});
    Dok              = Do(:,:,If{k});
    Dok              = distributed(Dok);
    Dok              = Dok(:);
    fh               = @(m) misfit_uq_slowness_smp(m,Dok,Q,modelk,paramsk);
    tic
    [h, Lpr, Lpm]    = fh(m);
    T(1) = toc;
    hh               = h;
    h                = double(h);
    Lpm              = diag(model.SIGMA_mp(:).^(-1));
    LAll             = {Lpm};
    betaAll          = {1};
    mmap             = velf(:);
    gammaAll         = {mmap(:)};
    hhh              = Hvs' * h * Hvs;
    C                = SmoothCov(model.n,a,b,c,dmax);
    Hfull            = hhh + inv(C);
    L                = chol(Hfull);
    %L =chol(hh);
    %figure;imagesc(reshape(1./diag(L),model.n));
end

Hvs              = opDiag(fun_df(vel(:)));

J                = oppJpenApp(hh.mt,hh.D,hh.Q,hh.model,hh.params,hh.G,hh.V,hh.EIG,hh.lambda,hh.AA,hh.Px,hh.distribute);

J2               = J * Hvs;
betaAll          = {1};
Lpm              = chol(C);
Lpm              = opInverse(Lpm');
JJ               = [J2; betaAll{1}*Lpm];

rm               = randn(size(JJ,1),nsmp);
tic
for i = 1:size(rm,2)
    dm(:,i) = lsqr(JJ,rm(:,i),1e-6,2000);
end
T(2) = toc;

for i = 1:nsmp
    SMP_RTO(:,i) = velf(:) + dm(:,i);
end

for i = 1:prod(model.n)
    std_RTO(i) = std(SMP_RTO(i,:));
end

mean_RTO = mean(SMP_RTO,2);
std_RTO  = reshape(std_RTO,model.n);
mean_RTO = reshape(mean_RTO,model.n);
SMP      = SMP_RTO;

WriteAllData('STD.mat',std_RTO,model.n,model.d,model.o);
WriteAllData('MEAN.mat',mean_RTO,model.n,model.d,model.o);
save('SMP.mat','SMP','-v7.3');



cd(curdir)
