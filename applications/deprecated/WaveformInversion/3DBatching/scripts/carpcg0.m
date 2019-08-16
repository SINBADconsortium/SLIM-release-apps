setpath;

expdir = [resultsdir '/' mfilename];
vfile  = [datadir '/overthrust.rsf']; 

p     = [16 8 4];
ns    = 1;
bs    = 1;
maxit = 2000;
tol   = 1e-6;

for k = 1:length(p)
    runCARPCG(expdir,vfile,p(k),ns,bs,tol,maxit);
end

