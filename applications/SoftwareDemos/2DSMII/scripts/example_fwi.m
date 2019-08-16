% This script runs the FWI example. For parallel runs; use 3 workers.
% The script may take several hours.

% directoy stuff
setpaths;
curdir = pwd;
expdir = [resultsdir '/fwi'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% model
[v,o,d,n]  = odnread([datadir '/bg_vp.odn']);
[v0,o,d,n] = odnread([datadir '/bg_v0.odn']);
mref       = 1e6./v.^2;
m0         = 1e6./v0.^2;

% model params
model.o    = o; model.d = d; model.n = n;
model.nb   = [20 20];
model.freq = 2.5:.5:20;
model.zsrc = 20;
model.xsrc = 0:100:7000;
model.zrec = 10;
model.xrec = 0:25:7000;
model.f0   = 15;
model.t0   = 0;

nsrc  = length(model.xsrc);
nrec  = length(model.xrec);
nfreq = length(model.freq);

Q = speye(nsrc);

% make data
D = F(mref,Q,model);
D = invvec(D,[nrec*nsrc nfreq]);

% inversion
If = {[1:3],[7:9],[19:21],[31:33]};

odnwrite('m_ls_0.odn',m0,o,d,n);

for k = 1:length(If)
    modelk      = model;
    modelk.freq = model.freq(If{k});
    Dk          = vec(D(:,If{k}));
    
    fh = @(x)f_ls(x,Q,Dk,modelk);
    
    options.fid     = fopen(['iter_ls_' num2str(k) '.log'],'w');
    options.itermax = 10;
    
    m0 = odnread(['m_ls_' num2str(k-1) '.odn']);
    
    mn = lbfgs_wWolfe(fh,m0,options);
    
    odnwrite(['m_ls_' num2str(k) '.odn'],mn,o,d,n);   
end


