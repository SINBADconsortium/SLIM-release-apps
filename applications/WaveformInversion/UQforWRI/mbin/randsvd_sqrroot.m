function [IHh DD]= randsvd_sqrroot(H,L,opts)
%% Hh = randsvd_sqrroot(H,L,opts)
%  function to create the spot operator of Hs = (H + L^{-1}L^{-t})
%  Input:
%  H
%  L
%  opts -
%  .rank - objective rank
%  .tol  - tolerence
%
%  Output:
%  IHh  - spot operator of the inverse square root of (H+LL')

tol  = 1e-6;
maxrank = size(H,1);

if isfield(opts,'tol')
    tol = opts.tol;
end

if isfield(opts, 'rank')
    Hrank = opts.rank;
else
    error('No rank of H is setted, Please set a rank for H to opts.rank')
end

if isfield(opts, 'rank')
    maxrank = opts.maxrank;
end

if isnumeric(L)
    L = opMatrix(L);
end

[U D V] = randsvd(L'*H*L,Hrank,tol,maxrank);
% [U D V] = svd(L'*H*L);

V       = U;
V       = opMatrix(V);

rD      = size(D,1);
rH      = size(H,1);



If      = opNull(rH);
DI      = opDiag((diag(D) + ones(rD,1)).^(-1/2) - ones(rD,1));
IHh     = L * (V*DI*V' + If);

DD = diag(D);

% If = eye(rH);
% DI = diag((diag(D) + ones(rD,1)).^(-1/2) - ones(rD,1));
% IHh     = L * (V*DI*V' + If);
