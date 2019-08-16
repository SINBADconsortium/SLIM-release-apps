function [x,res,niter] = CARPCG_OLD(R,idx,q,x,n,options)
% serial CARPCG_OLD algorithm for band-storage matrices
%
% use:
%   [x,res] = CARPCG_OLD(R,idx,q,x0,n,options)
%
% input:
%   {R,idx} - normalized matrix in band storage format, see mat2R
%   q       - right hand side
%   x0      - initial guess
%   n       - size of domain
%   options.w     - relaxation parameter, default = 1.5
%   options.tol   - CG tolerance, default = 1e-6
%   options.maxit - max. iterations, default = length(q)
%   options.ns    - number of double sweeps per iterations, default = 1
%   options.n_threads - number of CARP blocks to use for an additional
%   parallelization over threads. Each MATLAB worker will split up its block
%   into options.n_threads sub-blocks
%
% output:
%   x   - solution  
%   res - residual vector    
%
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

% parse options
w     = 1.5;
tol   = 1e-6;
maxit = length(q);
ns    = 1;

if isfield(options,'w')
    w = options.w;
end
if isfield(options,'tol')
    tol = options.tol;
end
if isfield(options,'maxit')
    maxit = options.maxit;
end
if isfield(options,'ns')
    ns = options.ns;
end

%%% Note: the following section also appears in CARPCG.m in this directory.

% If supplied in the call to the function, that value of n_threads takes
% precedence.
if isfield(options, 'n_threads')
	n_threads = options.n_threads;
% Otherwise use the environment variable.
else
	n_threads = str2num(getenv('OMP_NUM_THREADS'));
end
% If after all that n_threads is still not a valid value, fall back to the
% default.
if not(n_threads >= 1)
	n_threads = 1;
end
%%% End of copy-and-pasted section.

% Make sure diagonals of matrix are stored as rows
if size(R,1)==length(idx)
	% do nothing
elseif size(R,2)==length(idx)
	R = R.';
else
	error('Dimensions of R and idx do not match.');
end

% construct r.h.s. and operator
b = dsweep(R,idx,w,0*x,q,ns,n_threads);
Afun = @(x)(x - dsweep(R,idx,w,x,0*b,ns,n_threads));

% start CG 
tol  = (tol*norm(b))^2;
res  = zeros(maxit,1);
k    = 0;

r = b - Afun(x);
p = r;
normr = norm(r)^2;
res   = normr;

while (k<=maxit)&&(normr>tol)    

    q = Afun(p);
    
%      alpha = normr/(p'*q); I think that the sweeps are bugged and produce the 
%   complex conjugate of the correct solution. Check this asap!  - Lago
    alpha = normr/(q'*p);

    
    x = x + alpha*p;
    r = r - alpha*q;
    
    normr = norm(r)^2;
    beta  = normr/res(end);
    
    p = r + beta*p;
    
    k = k + 1;
    res(k) = normr;
end

res   = sqrt(res)/norm(b);
niter = k;
end


function x = dsweep(R,idx,w,x,b,ns,n_threads)
for k=1:ns
	% forward sweep
	x = sweepR_mex(R,idx,x,b,w,1,n_threads);
	% backward sweep
	x = sweepR_mex(R,idx,x,b,w,-1,n_threads);
end
end
