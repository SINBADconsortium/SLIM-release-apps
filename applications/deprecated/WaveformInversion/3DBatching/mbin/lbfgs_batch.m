function [x] = lbfgs_batch(fh,x0,options)
% Simple L-BFGS method with Wolfe linesearch
%
% use:
%   [xn,info] = mylbfgs(fh,x0,options)
%
% input:
%   fh - function handle to misfit of the form [f,g] = fh(x,I,eps)
%        where I is an indexset of size min(b0 + k*beta,bmax) at iteration
%        k and eps is a tolerance alpha^k*eps0.
%        f is the function value, g is the gradient of the same size
%        as the input vector x. 
%   x0 - initial guess
%
%   options.itermax - max iterations [default 10]
%   options.tol     - tolerance on 2-norm of gradient [1e-6]
%   options.M       - history size [5]
%   options.fid     - file id for output [1]
%   options.write   - save iterates to disk [0]
% 
%
%   options.eps0    - initial tolerance for misfit [1e-4]
%   options.epsmin  - smallest tolerance [1e-4]
%   options.alpha   - controls decrease of tolerance max(epsk = alpha^k*eps0,epsmin) [1]
%   options.b0      - initial batchsize [bmax]
%   options.bmax    - max. batch size   
%   options.beta    - controls growth of batchsize bk = min(b0 + k*beta,bmax) [0]
%   options.seed    - random seed for picking indexset I [1]
%   options.redraw  - redraw indexset at each iteration. [0]
%
% output:
%   xn - final estimate
%
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

if nargin<3
    options = [];
end

% various parameters
if isfield(options,'bmax')
   bmax = options.bmax;
else
    error('options.bmax must be provided.');
end

itermax    = getoption(options,'itermax',5);
M          = getoption(options,'M',5);
fid        = getoption(options,'fid',1);
tol        = getoption(options,'tol',1e-6);
write      = getoption(options,'write',0);
linesearch = getoption(options,'linesearch','wolfe');
eps0       = getoption(options,'eps0',1e-4);
epsmin     = getoption(options,'epsmin',1e-4);
alpha      = getoption(options,'alpha',0);
b0         = getoption(options,'b0',bmax);
beta       = getoption(options,'beta',0);
seed       = getoption(options,'seed',1);
redraw     = getoption(options,'redraw',0);
    
% initialization
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);
          
n = length(x0);
converged = 0;
iter = 0;
x = x0;
S = zeros(n,0);
Y = zeros(n,0);

epsk = eps0;
bk   = b0;

% initial evaluation
Ik     = randperm(bmax);
Ik     = Ik(1:floor(bk));
[f,g]  = fh(x,Ik,epsk);

nfeval = bk;
fprintf(fid,'# iter, # eval, stepsize, f(x)       , ||g(x)||_2\n');
fprintf(fid,'%6d, %6d, %1.2e, %1.5e, %1.5e\n',iter,nfeval,1,f,norm(g));
if write
    dlmwrite(['x_' num2str(iter) '.dat'],x);
end
% main loop
while ~converged
    
    % compute search direction
    s = B(-g,S,Y);
    p = -(s'*g)/(g'*g);
    
    if (p < 0)
        fprintf(fid,'Loss of descent, reset history\n');
        S = zeros(n,0);
        Y = zeros(n,0);
        s = B(-g,S,Y);
    end
    
    % linesearch
    switch linesearch
        case 'wolfe'
            [ft,gt,lambda,lsiter] = wWolfeLS(@(x)fh(x,Ik,epsk),x,f,g,s);
        otherwise
            [ft,gt,lambda,lsiter] = ArmijoLS(@(x)fh(x,Ik,epsk),x,f,g,s);
    end
    nfeval = nfeval + lsiter*bk;
    
    % update
    xt = x + lambda*s;

    S = [S xt - x];
    Y = [Y gt - g];

    if size(S,2)>M
        S = S(:,end-M+1:end);
        Y = Y(:,end-M+1:end);
    end
    f = ft;
    g = gt;
    x = xt;
    
    % update bk and epsk
    epsk = max(alpha*epsk,epsmin);
    bk   = min(bk + beta,bmax);
    if redraw
        Ik   = randperm(bmax);
    end
    Ik    = Ik(1:floor(bk));  
    
    % new evaluation
    [f,g] = fh(x,Ik,epsk);
    iter  = iter + 1;
    
    fprintf(fid,'%6d, %6d, %1.2e, %1.5e, %1.5e\n',iter,nfeval,lambda,f,norm(g));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],x);
    end
    
    % check convergence
    converged = (iter>itermax)||(norm(g)<tol)||(lambda<tol);
    
end

end

function z = B(x,S,Y)
% apply lbfgs inverse Hessian to vector
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   z = B(x,S,Y,b0,Binit)
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
if M>0 
    a = (Y(:,end)'*S(:,end)/(Y(:,end)'*Y(:,end)));
else
    a = 1/norm(x,2);
end
z = a*q;

% second recursion
for k = 1:M
    beta = rho(k)*Y(:,k)'*z;
    z    = z + (alpha(k) - beta)*S(:,k);
end
end

function [ft,gt,lambda,lsiter] = ArmijoLS(fh,x0,f0,g0,s0)
    
lambda = 2;
c1 = 1e-2;
done = 0;
lsiter = 0;
while ~done
    if lsiter < 10
        lambda = lambda/2;
        fp = f0 + c1*lambda*g0'*s0;
        [ft,gt] = fh(x0 + lambda*s0);
        done = (ft<=fp); 
        lsiter = lsiter + 1;
    else
        lambda = 0;
        break;
    end
end
end

function [ft,gt,lambda,lsiter] = wWolfeLS(fh,x0,f0,g0,s0)
% Simple Wolfe linesearch, adapted from
% (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3).
%
%

lsiter = 0;
c1 = 1e-2;
c2 = 0.9;
done = 0;
mu = 0;
nu = inf;
lambda = .5;
while ~done
    if nu < inf
        lambda = (nu + mu)/2;
    else
        lambda = 2*lambda;
    end
    
    if lsiter < 10
        [ft,gt] = fh(x0 + lambda*s0);
        lsiter = lsiter + 1;
    else
        lambda = 0;
        break;
    end
    
    %fprintf(1,'      >%d, %1.5e, %1.5e, %1.5e\n',lsiter, lambda, ft, gt'*s0);
    
    if ft > f0 + c1*lambda*g0'*s0
        nu = lambda;
    elseif gt'*s0 < c2*g0'*s0
        mu = lambda;
    else
        done = 1;
    end
end
end
