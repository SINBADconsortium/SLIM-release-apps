%% 
% This script implements the Stochastic LBFGS-B method described in 
%
% 'A unified 2D/3D software environment for large scale time-harmonic full waveform inversion' Curt Da Silva and Felix Herrmann
%
% applied to a 2D BG Compass example.  
%
% Each frequency batch should run in under an hour. 
%
% See the README file for instructions on how to download pre-run results
% and data.
%
% If you want to run in parallel, use 4 workers. It is not recommended to run this program in serial.
% The results are stored in the path defined in the script setpath.m
%
% The results are displayed in <bg2_batch_results.html bg2_batch_results.m>


% set parameters for experiment
label    = 'bg2_batch';   % name
vfile    = 'bg2v.rsf';    % reference velocity model 
v0file   = 'bg2v0.rsf';   % initial velocity model
datafile = 'bg2data.rsf'; % input data

% parameters for modeling
npml     = 20;
zsrc     = 20;
zrec     = 10;
f0       = 10;
t0       = 0;

% indices of frequency bands, total 17 bands of 3 frequencies each.
nf = 35;
batch_size = 3;
overlap = 1;
If = partition(nf,batch_size,overlap);

% min and max offset to use.
hmin     = 100;
hmax     = 3000;

% parameters for optimization
maxiter  = 5; % max. outer iterations per frequency band
subprob  = 1; %scaling factor, so each outer iteration 
              %is this number times as expensive as a full sweep through the data
N0       = 10;% initial batchsize on each worker
seed     = 1; % random seed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual script, you should not need to change anything below

setpath;
expdir = [resultsdir label '/'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
curdir = pwd;

% read models
[v,n,d,o]  = rsf_read_all([datadir vfile]);
[v0,n,d,o] = rsf_read_all([datadir v0file]);
mref       = 1e6./v(:).^2;
m0         = 1e6./v0(:).^2;

% read data
[D,nd,dd,od] = rsf_read_all([datadir datafile]);
D            = reshape(D,nd(1),nd(2),nd(3));

% model params
model.o    = o; model.d = d; model.n = n;
model.freq = od(3) + (0:nd(3)-1)*dd(3);
model.zsrc = zsrc;
model.xsrc = od(2) + (0:nd(2)-1)*dd(2);
model.zrec = zrec;
model.xrec = od(1) + (0:nd(1)-1)*dd(1);
model.f0   = f0;
model.t0   = t0;
model.unit = 's2/km2';

Q = speye(nd(2));

% offset mask
hh = model.xrec'*ones(1,nd(2)) - ones(nd(1),1)*model.xsrc;
params = struct; params.wri = false;
params.batch_mode = true;
params = default_fwi_params2d(params);
C = ones(nd(1),nd(2));
C(abs(hh)<=hmin) = 0;
C(abs(hh)>=hmax) = 0;
params.pdefunopts.offset_mask = C;
params.pdefunopts.helm_pml_max = npml;
nsrc = length(model.xsrc); nrec = length(model.xrec); nfreq = length(model.freq);

minv = min(vec(mref)); maxv = max(vec(mref));

% Inversion
modelk = model;
for k = 1:length(If)
    Ifk = If(k,:);
    disp(['Processing frequencies ' num2str(model.freq(Ifk))]);
    % select frequency band
    nf = length(Ifk);
    Dk = D(:,:,Ifk); 
    Dk = reshape(Dk,nrec,nsrc*nf);
    spmd,
        Dk = codistributed(Dk);
        Dk = redistribute(Dk,codistributor1d(2,codistributor1d.unsetPartition,size(Dk)));
    end    
    modelk.freq = model.freq(Ifk);
    
    % function handle to misfit
    fh = misfit_setup(m0,Q,Dk,modelk,params);
       
    % run inversion
    frugal_opts = struct;
    spmd, bmax = size(getLocalPart(Dk),2); end
    frugal_opts.bmax = [bmax{:}];
    frugal_opts.b0 = N0*ones(size(frugal_opts.bmax));
    frugal_opts.maxIter = maxiter;           
    frugal_opts.innerit = subprob; 
    mn = minfunc_batch(fh,m0,minv,maxv,frugal_opts);
    
    % write results
    rsf_write_all([expdir 'mn_' num2str(k) '.rsf'],{'out=stdout'},reshape(mn,n),d,o);
    
    m0 = mn;
end