function [x] = lbfgsbatch(fh,x0,options)
% Simple L-BFGS method with batching
%
% use:
%   [xn,info] = lbfgsbatch(fh,x0,options)
%
% input:
%   fh - function handle to misfit of the form [f,g] = fh(x,I)
%        where f is the function value, g is the gradient of the same size
%        as the input vector x and I is an indexset
%   x0 - initial guess
%
%   options.itermax - max iterations (10)
%   options.evalmax - max number of function evaluations (lsitermax*itermax*m)
%   options.tol     - tolerance on 2-norm of gradient (1e-6)
%   options.fid     - file id for output (1)
%   options.write   - write iterates to disk (false)
%   options.M       - history size (5)
%   options.lsitermax - max linesearch backtracks (10)
%   options.{c,beta}  - Armijo linesearch parameters (1e-2,.5)
%
%   for the batching options, the function is evaluated on part of the
%   index fh(x,I). The length of I grows as mk = gamma*mk + delta
%
%   options.m       - max size of indexset. REQUIRED
%   options.m0      - initial size of indexset (1)
%   options.gamma   - controls growth of batch (1)
%   options.delta   - controls growth of batch (1)
%   options.redraw  - whether to reshuffle the indices at each iteration (1)
%   options.scale   - scale values by current batchsize (true)
%   options.seed    - random seed to initialize randstream (randn)
%
% output:
%   xn - final estimate
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% various parameters
c         = 0.01;
beta      = 0.5;
M         = 5;
fid       = 1;
itermax   = 10;
evalmax   = 100;
lsitermax = 20;
tol       = 1e-6;
write     = 0;
% batching options
gamma     = 1;
delta     = 1;
m0        = 1;
redraw    = 1;
seed      = rand;
scale     = 1;

if isfield(options,'itermax')
    itermax = options.itermax;
end
if isfield(options,'lsitermax');
    lsitermax = options.lsitermax;
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
if isfield(options,'c')
    c = options.c;
end
if isfield(options,'beta')
    beta = options.beta;
end
if isfield(options,'gamma')
    gamma = options.gamma;
end
if isfield(options,'delta')
    delta = options.delta;
end
if isfield(options,'m')
    m = options.m;
else
    fprintf(fid,'max batchsize options.m is a required parameter!\n');
    return;
end
if isfield(options,'m0')
    m0 = options.m0;
end
if isfield(options,'redraw')
    redraw = options.redraw;
end
if isfield(options,'seed')
    seed = options.seed;
end
if isfield(options,'scale')
    scale = options.scale;
end
evalmax = 2*itermax*m;
if isfield(options,'evalmax')
    evalmax = options.evalmax;
end

% initialization
n = length(x0);
converged = 0;
iter = 0;
x = x0;
S = zeros(n,M);
Y = zeros(n,M);
mk = m0;
nfeval = 0;
cnt = 0;

rndstrm = RandStream.create('mt19937ar','Seed',seed);
I  = randperm(rndstrm,m);

% main loop
while ~converged
    [f,g]  = fh(x,I(1:mk));
    nfeval = nfeval + mk;
    if scale
        f = m*f/mk;
        g = m*g/mk;
    end
    
    % compute search direction
    s = B(-g,S(:,max(M+1-iter,1):end),Y(:,max(M+1-iter,1):end));

    % Armijo linesearch
    lambda0 = 1;
    if iter==0
        lambda0 = 1/max(norm(g,1),1);
    end
    lambda = lambda0/beta;ft = 1; fp = 0;lsiter = 0;
    while ft > fp
        lambda = lambda*beta;
        lsiter = lsiter + 1;
        
        xt = x + lambda*s;

        [ft,gt] = fh(xt,I(1:floor(mk)));
        nfeval  = nfeval + floor(mk);
        if scale
            ft = m*ft/floor(mk);
            gt = m*gt/floor(mk);
        end

        fp = f + c*lambda*g'*s;
       
        if lsiter > lsitermax
            fprintf(fid,'Linesearch failure, quitting\n');
            return
        end
    end
    
    % print result
    fprintf(fid,'%5d, %5d, %5d, %10.5e, %10.5e, %10.5e\n',iter,nfeval,mk,lambda,ft,norm(gt));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],xt);
    end
    
    % update
    S = [S(:,2:end) xt - x];
    Y = [Y(:,2:end) gt - g];

    x = xt;

    % increase batchsize
    mk = min(gamma*mk + delta,m);
    if redraw
        I = randperm(rndstrm,m);
    end
    
    iter = iter + 1;
    
    % check convergence
    converged = (nfeval>=evalmax)||(iter>=itermax)||(norm(g)<tol)||(norm(S(:,end))<tol);
    
end

end

function z = B(x,S,Y)
% apply lbfgs inverse Hessian to vector
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   z = B(x,S,Y)
%
% input:
%   x - vector of length n
%   S - history of steps in n x M matrix
%   Y - history of gradient differences in n x M matrix
%
% output
%   z - vector of length n
%

M = size(S,2);

alpha = zeros(M,1);
rho   = zeros(M,1);
for k = 1:M
    rho(k) = 1/(Y(:,k)'*S(:,k));
end
q = x;
% first recursion
for k = M:-1:1
    alpha(k) = rho(k)*S(:,k)'*q;
    q        = q - alpha(k)*Y(:,k);
end

% apply `initial' Hessian
if M==0
    B0 = 1;
else
    B0 = Y(:,end)'*S(:,end)/(Y(:,end)'*Y(:,end));
end
z  = B0*q;

% second recursion
for k = 1:M
    beta = rho(k)*Y(:,k)'*z;
    z    = z + (alpha(k) - beta)*S(:,k);
end
end



