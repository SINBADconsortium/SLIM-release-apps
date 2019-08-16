setpath;

expdir = [resultsdir '/' mfilename];
vfile  = [datadir '/overthrust.rsf']; 

p     = [4];
ns    = 1;
bs    = 1;
maxit = 2000;
tol   = 1e-6;

runCARPCG(expdir,vfile,p,ns,bs,tol,maxit);

system(['cat ' expdir '/table.dat']);
