function runFWI(expdir, vfile, datafile, params)
% run 3D FWI with batching
%
% use
%   runFWI(expdir, vfile, datafile, params)
%
% input:
%   expdir   - directory for output
%   vfile    - .rsf file with initial model (s^2/km^2)
%   datafile - .rsf file with data
%   params.nb      - # of points for PML (default = 10)
%          hmin    - min. offset to use (default = 0)
%          maxiter - # of L-BFGS iterations (default = 10)
%          {eps0,epsmin,alpha} - control accuracy of PDE solves (epsk = max(alpha^k*eps0,epsmin)) (default = {1e-4,1e-4,0})
%          {b0,beta}           - control batchsize (bk = min(b0 + beta*k,nsrc)) (default = {inf,0})
%          seed                - seed for random generator (default = 1)
%          redraw              - 1: redraw batch, 0: don't redraw (default = 0)
%          If                  - cell array with frequency bands for multiscale inversion         

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

nb      = getoption(params,'nb',10);
hmin    = getoption(params,'hmin',0);
maxiter = getoption(params,'maxiter',10);
eps0    = getoption(params,'eps0',1e-4);
epsmin  = getoption(params,'epsmin',1e-4);
alpha   = getoption(params,'alpha',0);
b0      = getoption(params,'b0',inf);
beta    = getoption(params,'beta',0);
seed    = getoption(params,'seed',1);
redraw  = getoption(params,'redraw',0);
If      = getoption(params,'If',{1});
k0      = getoption(params,'k0',1);
%%
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% read model
[m0,n,d,o]      = rsf_read_all(vfile);

% read data
[data,nd,dd,od] = rsf_read_all(datafile);

% setup parameters
model.o  = o;
model.d  = d;
model.n  = n;
model.nb = nb*[1 1 1];

[model.zrec,model.xrec,model.yrec,model.zsrc,model.xsrc,model.ysrc,model.freq] = odn2grid(od,dd,nd);

model.f0   = 0;
model.t0   = 0;

% inversion
Q  = speye(prod(nd(4:6)));

options.srcest = 1;
options.hmin   = hmin;
options.maxit  = prod(n);

rsf_write_all('m0.rsf',{},reshape(m0,n),d,o);
modelk = model;
for k = k0:length(If)

    m0 = rsf_read_all(['m' num2str(k-1) '.rsf']);
    m0 = m0(:);
    
    modelk.freq = model.freq(If{k});
    
    fh = @(x,I,eps)misfit(x,Q,I,eps,vec(data(:,:,:,:,:,:,If{k})),modelk,options);

    opts.itermax = maxiter;
    opts.fid     = fopen(['iter' num2str(k) '.log'],'w');
    opts.write   = 1;
    opts.eps0    = eps0;
    opts.epsmin  = epsmin;
    opts.alpha   = alpha;
    opts.b0      = min(b0,prod(nd(4:6)));
    opts.bmax    = prod(nd(4:6));
    opts.beta    = beta;
    opts.seed    = seed;
    opts.redraw  = redraw;

    mn = lbfgs_batch(fh,m0,opts);

    rsf_write_all(['m' num2str(k) '.rsf'],{},reshape(mn,n),d,o);

    mkdir(['xm' num2str(k)]);
    movefile('*.dat',['./xm' num2str(k)]);
end
