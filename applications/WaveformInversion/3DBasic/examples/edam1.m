%% This is an example script for inverting the simple 3D Edam model 
% with transmission data. This is primarily to demonstrate how to
% use the new 3D forward modelling + FWI code with model coarsening
%
% Author: Curt Da Silva, 2015
%
cur_dir = pwd; base_dir = cur_dir(1:end-length('examples'));
results_dir = [base_dir 'results/'];
if exist(results_dir,'dir')==0, mkdir(results_dir); end

% Free parameters
% Number of points in the x,y,z directions
nx = 50; ny = 50; nz = 50;

% Background velocity
vlow        = 2000;

% Perturbation velocity
vhi         = 2200;

% Number of right hand sides over which to solve Helmholtz in
% parallel using multi-threading
numsrc_mt = 1;

% Number of parallel workers
num_workers = 4;

% Number of sources (for x, y dimensions individually)
% Total number of sources is this number^2
numsrc = 4;

% Number of LBFGS iterations per frequency batch
num_iter = 10;

% If true, coarsen model depending on max wavelength of each freq batch
coarsen_model = false;

% Linear solve tolerance for FWI
linearsolve_tol = 1e-6;

%% Model setup 
model       = struct;
x = linspace(0,3000,nx); y = linspace(0,3000,ny); z = linspace(0,3000,nz);
[o,d,n] = grid2odn(x,y,z);

model.o     = o;
model.d     = d;
model.n     = n;
model.unit  = 'm/s';

% Perturbation
[X,Y,Z] = ndgrid(1:length(x),1:length(y),1:length(z));
R = nx/4;
v = vlow*ones(model.n);
v( (X-nx/2).^2 + (Y-ny/2).^2 + (Z-ny/2).^2 < R^2 ) = vhi;

v0 = vlow*ones(model.n); v0(1) = 0.999999*v0(1);
if strcmp(model.unit,'s2/m2'), v0 = v0.^(-2); end

%% Set up source/receiver grid
% options for solving the helmholtz equations
lsopts = LinSolveOpts();
lsopts.tol = 1e-6;
lsopts.maxit = 2000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

% pdefunc options
pdeopts = PDEopts();
pdeopts.helm_dt = model.d;
pdeopts.debug_mode = false;
pdeopts.helm_pml_max = 10;
pdeopts.numcompsrc = numsrc_mt;

opts = struct;
opts.pdefunopts = pdeopts;
opts.lsopts = lsopts;
opts.subsample_model = coarsen_model;
opts.disp_progress = false;
opts.debug_mode = false;

[x,y,z] = odn2grid(model.o,model.d,model.n);

% Transmission experiment
model.xsrc = x(round(linspace(3,nx-3,numsrc)));
model.ysrc = y(round(linspace(3,ny-3,numsrc)));
model.zsrc = min(z)+model.d(3);
model.xrec = x(3:1:end-2);
model.yrec = y(3:1:end-2);
model.zrec = max(z)-model.d(3);

% source wavelet time shift
model.t0   = 0; 

% ricker wavelet peak frequency
model.f0 = 10;

% frequencies in Hz
model.freq = 2:4;

nrx = length(model.xrec); nry = length(model.yrec);  
nfreq = length(model.freq);
nsx = length(model.xsrc); nsy = length(model.ysrc);
nsrc = nsx*nsy;
nrec = nrx*nry;

Q = speye(nsrc);   
sz_data = [nrx,nry,nsx,nsy,nfreq];

%% Generate data + optimize
minfunc_opts = struct;
minfunc_opts.maxIter = num_iter;
minfunc_opts.optTol = -1;
minfunc_opts.obj_aux = false;

if parpool_size()==0, parpool_open(num_workers); end

Dobs = F3d(v,Q,model,opts);

%% Frequency continuation FWI + model coarsening

% user-specified tolerance for helmholtz solves
lsopts.tol = linearsolve_tol;

% work on a single frequency at a time
freq_batch = num2cell(1:nfreq,[1,nfreq]);

vhist = cell(length(freq_batch),1);
vest = vec(v0);
for j=1:length(freq_batch)
    % Select only sources at this frequency batch
    srcfreqmask = false(nsrc,nfreq);
    srcfreqmask(:,freq_batch{j}) = true;
    opts.srcfreqmask = srcfreqmask;
    [obj,vest_sub,comp_grid] = misfit_setup(vest,Q,Dobs,model,opts);    
    
    % Optimization on coarser grid
    vest_sub = minConf_TMP(obj,vest_sub,vlow,vhi,minfunc_opts);
    vest = comp_grid.to_fine*vest_sub;
    vhist{j} = vest;        
end

%% Save results + data residuals for visualization

save([results_dir 'edam_fwi.mat'],'v','v0','vest','model','vhist');

% Generate data for a single source/freq pair, just for plotting purposes
sx = 2; sy = 1; f = 1; 
Isrc = sub2ind([nsx,nsy,nfreq],sx,sy,f);
srcfreqmask = false(nsrc,nfreq);

srcfreqmask(Isrc) = true;
opts.srcfreqmask = srcfreqmask;

Dinit = F3d(v0,Q,model,opts);

D0 = gather(Dobs(:,Isrc));

Dfinal = F3d(vest,Q,model,opts);

init_res = reshape(gather(D0-Dinit),nrx,nry);
final_res = reshape(gather(D0-Dfinal),nrx,nry);

save([results_dir 'edam_fwi.mat'],'-append','init_res','final_res','sx','sy','f');
