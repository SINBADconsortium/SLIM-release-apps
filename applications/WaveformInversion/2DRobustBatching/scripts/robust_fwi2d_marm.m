%%
% The |mbase|, |mlsls|, |mstls| and |mstst| scripts reproduce the results
% from A.Y. Aravkin, T. van Leeuwen and F.J. Herrmann - Source estimation
% for frequency-domain FWI with robust penalties, EAGE Expanded abstracts
% 2012, submitted.
%
% It is not advisable to run these in serial and/or interactive mode since
% it will take several hours to run. 
%
% See the README file for intstructions on how to download pre-run results
% and data.
%
% If you want to run in parallel, use a divisor of 12 workers. 
% The results are stored in the path defined in the script setpath.m
% 
%
% robust_fwi2d_marm(p,fh1,fh2,label);
%   p - percentage of outliers
%   fh1  - function handle to compute residual
%   fh2  - function handle to compute source estimation
%   label - experiment label 

function robust_fwi2d_marm(p,fh1,fh2,label)
    
% set parameters for experiment
vfile    = 'marmv.rsf';    % reference model
v0file   = 'marmv0.rsf';   % initial model
datafile = 'marmdata.rsf'; % input data

% modeling parameters
nb       = 50;
zsrc     = 10;
zrec     = 10;
f0       = 10;
t0       = 0;

% parameters for generating outliers in the data 
a        = 1e6;  % amplitude of noise
seed     = 1;    % random see

% function handles to misfits used for residual and source estimation


% number of iterations for L-BFGS
maxiter  = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual script, you should not need to change anything below

% directory stuff
setpath;
expdir = [resultsdir label '/'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end

% model
[v,n,d,o]  = rsf_read_all([datadir vfile]);
[v0,n,d,o] = rsf_read_all([datadir v0file]);
mref       = 1e6./v(:).^2;
m0         = 1e6./v0(:).^2;

% read data
[D,nd,dd,od] = rsf_read_all([datadir datafile]);
D            = reshape(D,nd(1)*nd(2),nd(3));

% model params
model.o    = o; model.d = d; model.n = n;
model.nb   = nb.*[1 1];
model.freq = od(3) + [0:nd(3)-1]*dd(3);
model.zsrc = zsrc;
model.xsrc = od(2) + [0:nd(2)-1]*dd(2);
model.zrec = zrec;
model.xrec = od(1) + [0:nd(1)-1]*dd(1);
model.f0   = f0;
model.t0   = t0;

Q = speye(nd(2));

% generate outliers by replacing a fraction p of the data with
% random normal (0, a)
s = RandStream.create('mt19937ar','seed',seed);
I = randi(s,prod(nd),floor(p*prod(nd)),1);
D(I) = a*randn(s,size(I));

D = distributed(vec(D));

% inversion
params.misfit  = fh1;
params.src_est = fh2;

fh = @(x)JI_src(x,Q,D,1:nd(2),model,params);

options.maxIter = maxiter;

mn = minConf_TMP(fh,m0,min(mref),max(mref),options);

% final evaluation
[fn,gn,wn] = fh(mn);

% write results
rsf_write_all([expdir 'wn.rsf'],{'out=stdout'},wn,[dd(3) 0],[od(3) 0]);
rsf_write_all([expdir 'mn.rsf'],{'out=stdout'},reshape(mn,n),d,o);

end