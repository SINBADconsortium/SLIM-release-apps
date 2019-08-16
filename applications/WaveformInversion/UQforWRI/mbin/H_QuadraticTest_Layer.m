% resultsdir = ['./test'];
% mkdir(resultsdir)
% cd(resultsdir);
addpath(genpath('/home/slim/zfang/Project/SofterRelease/SLIM-release-apps/tools/functions'));
addpath(genpath('/home/slim/zfang/Project/SofterRelease/SLIM-release-apps/tools/algorithms'));
addpath('/home/slim/zfang/Project/SofterRelease/SLIM-release-apps/tools/solvers/Krylov');
addpath('/home/slim/zfang/Project/Paper/Journal/SLIM.Papers.UQofWRI/code');
%vfile  = '/scratch/slim/zfang/Result/WRI_SRCEST/BGModel/TrueSRC_20HzCent/m10.mat';
%vtfile = '/scratch/slim/zfang/Model/BG2D/BG450kmTrue.mat';
vfile  = '/scratch/slim/zfang/Result/UQWRI/LayerModel/Paper/NoiseTest/m10.mat';
%vfile = '/scratch/slim/zfang/Model/Layer/Layer1_true_3km.rsf';
vifile = '/scratch/slim/zfang/Model/Layer/Layer1_ini_31.mat';
vtfile = '/scratch/slim/zfang/Model/Layer/Layer1_true_3km.rsf';
ResDir = '/scratch/slim/zfang/Result/UQWRI/LayerModel/Paper/QuadraticApp';

if ~exist(ResDir,'dir')
    mkdir(ResDir);
end

copyfile('./H_QuadraticTest_Layer.m',ResDir);

Mstr   = 'Vn';
lexp   = 0;



niter = 1000;
tol      = 1e-6;
dmratio  = 1;

% model
[vi n d o] = ReadAllData(vifile);
[vt]       = ReadAllData(vtfile);
[z x]      = odn2grid(o,d,n);

% acquisition
nfreq = 2;
fmin  = 21;fmax = 22;

model.o = o;
model.d = d;
model.n = n;
model.nb = [20 20];
%model.xrec = [10:10:1980];
model.xrec = [1000];
model.zrec = 20;
%model.xsrc = [80:80:1980];
model.xsrc = [1000];
model.zsrc = 20;

model.f0   = 0;
model.t0   = 0;
model.freq = linspace(fmin,fmax,nfreq);
nsrc = length(model.xsrc);

% data
m = 1e0./vt(:).^2;
Q = speye(nsrc);
D = F(m,Q,model);

% inversion
mi    = 1./vi(:).^2;
m0    = m;


params.C          = 1;
params.srcest     = 0;
params.lambda     = 10^lexp;
options.itermax   = 10;
params.wri        = 1;
params.dist_mode  = 'freq';
params.hessian    = 'hess-app';%'sparse';
params.Threads    = 10;
params.distribute = 1;
params.Hvel       = 0;

    
%fh_red = @(x)misfit_red(x,Q,D,model,params);
%fh_pen = @(x)misfit_pen(x,Q,D,model,params);

fh = misfit_setup(Q,D,model,params);

[f g H] = fh(m0);
disp('start to compute Diag');
%Hdiag = Diags(H);

dm    = randn(size(m0));
dm    = dm / norm(dm) * norm(m0) * dmratio;
%dm    = m0 - mi;

h     = [-1:0.2:1];
%g     = zeros(size(dm));
funq  = @(x) f+g'*(x)+0.5*x'*(H*x);
%funqd = @(x) f+g'*(x)+0.5*x'*(Hdiag.*x);

for i = 1:length(h)
    mi      = m0 + h(i) * dm;
    ft(i)   = fh(mi);
    fq(i)   = funq(h(i)*dm);
%    fqd(i)  = funqd(h(i)*dm);
end

%figure;plot([ft(:) fq(:) fqd(:)]);
figure;plot([ft(:) fq(:)]);
%legend('Ture misfit', 'Quadratic apprximation with approximated Hessian', 'Quadratic approximation with diagonal of the Hessian');
legend('Ture misfit', 'Quadratic apprximation with approximated Hessian');

Result.ft   = ft;
Result.fq   = fq;
Result.fqd  = fqd;
Result.h    = h;
Result.dmr  = dmratio;

filestr     = [ResDir '/' 'Quadratic_' Mstr '_' num2str(lexp) '.mat'];
save(filestr,'Result');



%figure;plot([ft(:) fq(:)]);

% hh  = 3.^[-3:-1:-8];
% funqd2 = @(x) f+g'*x+0.5*x'*(H*x);
% for i = 1:length(hh)
%     mi       = m0 + hh(i) * dm;
%     ft2(i)   = fh(mi);
%     fq2(i)   = funqd2(hh(i)*dm);
%     df(i)    = abs(ft2(i) - fq2(i));
% end
% figure;loglog(hh,df);hold on; loglog(hh,hh.^3);
%
%
%
%
% keyboard
% HH = h + opNull(size(h,1)) * 0.1;
%
% IH = oppIHpenApp(HH, niter, tol);
%
% keyboard
%
% %mr = lbfgs_wWolfe(fh_red,m0,options);
%
% %mp = lbfgs_wWolfe(fh_pen,m0,options);
%
% mlbfgs  = lbfgs_wWolfe(fh,m0,options);
% opt.solver = 'backslash';
% opt.maxiter = 10;
% mnewton = newton(fh,m0,1e-6,opt);
%
%
% % plot results
% vlbfgs  = reshape(1e3./sqrt(mlbfgs),n);
% vnewton = reshape(1e3./sqrt(mnewton),n);
%
% figure;
% subplot(1,2,1);
% imagesc(x,z,vlbfgs);
% a = caxis;
% subplot(1,2,2);
% imagesc(x,z,vnewton)
% caxis(a)
