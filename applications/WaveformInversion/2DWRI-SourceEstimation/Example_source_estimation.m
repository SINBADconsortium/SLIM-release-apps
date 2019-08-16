
%% Wavefield Reconstruction Inversion with Source Estimation
%This is a synthetic example of WRI with source estimation. We mimic the
%transmission case. Sources are loacted at the left side of the model while
%receivers are located at the right side of the model.
%In this example, we assume that the source functions are unknown and we will
%estimate them on-the-fly.

% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: March, 2018.

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (zfang@eoas.ubc.ca)

%This script will generate data and invert using the same modeling code.
%Model is the Camembert.

%% startup
close all
startup

%results directory
basedir    = pwd;
resultsdir = [basedir,'/results/'];

if ~isdir(resultsdir)
    mkdir(resultsdir);
end
cd(resultsdir)

%% Set up the model.

n = [101, 101];
d = [20, 20];
o = [0, 0];

vel    = 2*ones(n);
velini = vel;
for i = 1:n(1)
  for j = 1:n(2)
    if norm([i-51,j-51])<18
        vel(i,j) = 2.25;
    end
  end
end

figure;imagesc(vel)

model.o    = o;
model.d    = d;
model.n    = n;
model.f0   = 6;
model.t0   = 0;
model.zsrc = d(1):d(1)*2:d(1)*(n(1)-2);
model.xsrc = d(1)*2;
model.zrec = d(1):d(1):d(1)*(n(1)-2);
model.xrec = d(1)*(n(1)-2);
model.freq = 2:4:10;
model.nb   = [20 20];

params.wri       = 1;
params.hessian   = 'sparse';
params.nthreads  = 6;
params.dist_mode = 'freq';
params.lambda    = 1;
params.extend    = 1;
params.SrcEst    = 0;

if isfield(params,'pdefunopts')==0
    params.pdefunopts = PDEopts2D();
end
if isfield(params,'lsopts')==0
    solve_opts = LinSolveOpts();
    solve_opts.solver = LinSolveOpts.SOLVE_LU;
    params.lsopts = solve_opts;
end

options.itermax = 1000;
options.M       = 20;
options.write   = 0;

mt = vec(1./vel.^2);
mi = vec(1./velini.^2);

%% Set up the generation of synthetic data in the true velocity model
model.nrec = length(model.zrec);
model.nsrc = length(model.zsrc);

% source matrix
Q = speye(model.nsrc);

% generate data
fprintf('Start to generate the data\n');
D           = F(mt,Q,model,params);
Ds          = gather(D);
fprintf('Data is generated\n');

fh          = misfit_setup(Q,D,model,params);

fprintf('Start the inversion with the true source function\n');
[mf_t info_t] = lbfgs_wWolfe(fh,mi,options);
vf_t          = reshape(1./sqrt(mf_t),model.n);
fprintf('The inversion with the true source function finishes\n');
figure;imagesc(reshape(vel,model.n));caxis_a=caxis;
figure;imagesc(reshape(vf_t,model.n));caxis(caxis_a);

Result.mf   = mf_t;
Result.o    = o;
Result.d    = d;
Result.n    = n;
Result.info = info_t;
save('./Result_truesrc.mat','Result');

model_old        = model;

params.SrcEst    = 1;
model.f0         = 0;
fh               = misfit_setup(Q,D,model,params);
fprintf('Start the inversion with the true source function\n');
[mf_e info_e]    = lbfgs_wWolfe(fh,mi,options);
vf_e             = reshape(1./sqrt(mf_e),model.n);
fprintf('The inversion with the true source function finishes\n');
figure;imagesc(reshape(vf_e,model.n));caxis(caxis_a);
Result.mf   = mf_e;
Result.o    = o;
Result.d    = d;
Result.n    = n;
Result.info = info_e;
save('./Result_estsrc.mat','Result');

[f g h w faux] = fh(mf_e);
D1             = F(mf_e,Q,model,params);
D              = gather(D);
D1             = gather(D1);
D              = reshape(D,model.nrec,model.nsrc,length(model.freq));
D1             = reshape(D1,model.nrec,model.nsrc,length(model.freq));
for i = 1:length(model.freq)
  D1(:,:,i) = D1(:,:,i) * faux.est_src(1,i);
end

save('./Dobs.mat','D')
save('./Dpred.mat','D1')

src_est  = faux.est_src(25,:);
src_true = fwi_wavelet(model.freq,model_old.t0,model_old.f0);

save('./src_est.mat','src_est');
save('./src_true.mat','src_true');



cd(basedir)
