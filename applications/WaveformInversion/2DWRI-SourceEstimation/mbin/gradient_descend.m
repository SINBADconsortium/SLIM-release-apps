function x = gradient_descend(fh,x,options)

% various parameters
M         = 5;
fid       = 1;
itermax   = 10;
tol       = 1e-6;
write     = 0;
init_step = 0;

if isfield(options,'itermax')
    itermax = options.itermax;
end
if isfield(options,'M');
    M = options.M;
end
if isfield(options,'fid')
    fid = options.fid;
end
if isfield(options,'write')
    write = options.write;
end
if isfield(options,'tol')
    tol = options.tol;
end
if isfield(options,'stepsize')
    stepsize = options.stepsize;
end

[f,g]  = fh(x);
nfeval = 1;
iter   = 0;
fprintf(fid,'# iter, # eval, stepsize, f(x)       , ||g(x)||_2\n');
fprintf(fid,'%6d, %6d, %1.2e, %1.5e, %1.5e\n',iter,nfeval,1,f,norm(g));
if write
    dlmwrite(['x_' num2str(iter) '.dat'],x);
end
iter = 1;

for i = 1:itermax
    x = x - stepsize * g;
    [f,g]  = fh(x);
    fprintf(fid,'%6d, %6d, %1.2e, %1.5e, %1.5e\n',iter,nfeval,1,f,norm(g));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],x);
    end
    iter = iter+1;
end




