function d = TraceNorm_dual(x,weights, params)

% dual of trace norm is operator norm i.e maximum singular value

E = reshape(x,params.numr,params.numc);

d = svds(E,1);

