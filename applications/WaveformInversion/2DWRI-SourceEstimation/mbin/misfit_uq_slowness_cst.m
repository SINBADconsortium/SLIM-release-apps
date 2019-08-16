function [f g] = misfit_uq_slowness_cst(m,Dok,Q,modelk,paramsk,nlayercst,mcst)
m = reshape(m,modelk.n(1)-nlayercst,modelk.n(2));
m = [ones(nlayercst,modelk.n(2))*mcst; m];
m  = m(:);
[f g] = misfit_uq_slowness(m,Dok,Q,modelk,paramsk);
g       = reshape(g,modelk.n);
g       = g(nlayercst+1:end,:);
g       = g(:);


