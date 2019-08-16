setpath;

expdir = [resultsdir '/' mfilename];
vfile  = [datadir '/overthrust.rsf']; 

p     = [8];
ns    = 50;
bs    = [1 2 5 50];
maxit = 1000;
tol   = 1e-6;

for k = 1:length(bs)
    runCARPCG(expdir,vfile,p,ns,bs(k),tol,maxit);
end

system(['cat ' expdir '/table.dat']);
