function [ x,hst ] = BICR( A,b,x,opts )
%BICR BICR algorithm for solving nonsymmetric linear systems
% from 'An extension of the conjugate residual method to nonsymmetric
% linear systems' - Sogabe, Sugihara, Zhang, 2007. Supports
% multiple right hand sides.
%
% Curt Da Silva, 2015
%
% Usage:
%
%   [x,hst] = BICR(A,b,x,opts);
%
% Input:
%    A      - matrix or spot operator 
%    b      - right hand side
%    x      - initial guess
%    opts   - options struct
%        .maxit       - maximum number of iterations (default: size(b,1));
%        .precond     - preconditioner (either a spot operator or a
%                       function handle, default: @(x) x);
%        .precond_adj - adjoint precondtiioner, if opts.precond is
%                       a function handle (default: @(x) x)
%        .tol         - relative residual tolerance (default: 1e-6)
% 

maxit       = check_field(opts,'maxit',size(b,1));
P           = check_field(opts,'precond',opDirac(size(A,1)));
tol         = check_field(opts,'tol',1e-6);
out_freq    = check_field(opts,'out_freq',0);
if isa(P,'function_handle')    
    precond_adj = check_field(opts,'precond_adj',@(x) x);
    P = opFunction(size(A,2),size(A,2),{P,precond_adj},1,1);
end
A = A*P;
nrhs = size(b,2);

if ~isempty(x) && nnz(x)>0
   r = b-A*x;
else
   r = b; 
end
rs = r;
Ar = A*r; Ars = A'*rs;
p = zeros(size(r)); ps = zeros(size(r));
Ap = 0*p; Aps = 0*ps;
beta = spdiags(zeros(nrhs,1),0,nrhs,nrhs);
itr = 1;
if nrhs==1
    innerprod = @(x,y) vec(x)'*vec(y);
else
    innerprod = @(x,y) sum(conj(x).*y,1);
end

hst = zeros(maxit,1); hst(1) = max(norm(r));
while true                   
    p = r + p*beta;
    ps = rs + ps*conj(beta);
   
    Ap = Ar + Ap*beta; Aps = Ars + Aps*conj(beta);
    
    rAr = innerprod(rs,Ar);
    
    if nrhs==1
        alpha = rAr/innerprod(Aps,Ap);
    else
        alpha = spdiags(vec(rAr./innerprod(Aps,Ap)),0,nrhs,nrhs);
    end
    x = x + p*alpha;
    
    r = r - Ap*alpha;
    rs = rs - Aps*conj(alpha);
    
    Ar = A*r; Ars = A'*rs;
    
    if nrhs==1
        beta = innerprod(rs,Ar)/rAr;
    else
        beta = spdiags(vec(innerprod(rs,Ar)./rAr),0,nrhs,nrhs);
    end
    
    itr = itr+1;    
    hst(itr) = max(norm(r));
    if out_freq > 0 && mod(itr-1,out_freq)==0
        
        disp(['k : ' num2str(itr,'%3d') ' res ' num2str(max(norm(r)./norm(b)),'%3.3e')]);
    end
    if itr >= maxit, break; end
    if max(norm(r)./norm(b)) < tol, break; end
end
hst = hst(1:itr);

end

