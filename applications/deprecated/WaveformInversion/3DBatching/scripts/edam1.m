setpath;

expdir   = [resultsdir '/' mfilename];
vfile    = [datadir '/edam_m0.rsf']; 
datafile = [datadir '/edam_data.rsf']; 

params.nb      = 10;
params.hmin    = 0;
params.maxiter = 10;
params.eps0    = 1e-3;
params.epsmin  = 1e-3;
params.alpha   = 1;
params.b0      = inf;
params.beta    = 0;
params.seed    = 1;
params.redraw  = 0;
params.If      = {[1:3]};

runFWI(expdir,vfile,datafile,params);