function [x] = BICGSTAB(A,b,x,opts)
% BICGSTAB algorithm with flexible preconditioning, implementation
% of Alg. 1 in 'ANALYSIS AND PRACTICAL USE OF FLEXIBLE BICGSTAB',
% Chen, Mcinnes, Zhang
%
% Curt Da Silva, 2015
% 
% Usage:
%   x = BICGSTAB(A,b,x0,opts);
% 
% Input:
%   A    - matrix or spot operator
%   b    - right hand side
%   x0   - initial guess (default: 0)
%   opts.out_freq - display output every # of iterations (default: 0)
%   opts.tol      - relative residual tolerance
%   opts.maxit    - max number of iterations
%   opts.precond  - either a function handle (preconditioner for the
%                  forward problem) or a SPOT operator (both forward + adjoint prec)
%
%       .precond_adj - if opts.precond is a function handle, the
%                      preconditioner for the adjoint problem
%   
% Output:
%   x    - solution to Ax = b, within the specified tolerance
%
    maxit       = check_field(opts,'maxit',size(b,1));
    precond     = check_field(opts,'precond',@(x)x);
    tol         = check_field(opts,'tol',1e-6);
    out_freq    = check_field(opts,'output_freq',0);
    if isa(precond,'opSpot')
        P = @(x) precond*x;
    else
        P = precond;
    end

    if ~isempty(x) && nnz(x)>0
        r = b-A*x;
    else
        r = b; 
        x = zeros(size(b));
    end
    r0_bar = r;
    normb = norm(b);
    p = r;
    itr = 0;
    
    innerprod = @(x,y) x'*y;    
    while true
        p_tilde = P(p);
        Ap_tilde = A*p_tilde;
        alpha = innerprod(r,r0_bar) ./innerprod(Ap_tilde,r0_bar);
        s = r - Ap_tilde*alpha;
        s_tilde = P(s);
        As_tilde = A*s_tilde;
        omega = innerprod(As_tilde,s)./innerprod(As_tilde,As_tilde);
        xnew = x + p_tilde*alpha + s_tilde*omega;
        rnew = s - As_tilde*omega;
        beta = innerprod(rnew,r0_bar)./innerprod(rnew,r0_bar) .* alpha ./ omega;
        pnew = rnew + (p - Ap_tilde*omega) * beta;
        
        p = pnew; r = rnew; x = xnew;
        if max(norm(r)./normb) < tol
            break;
        end         
        
        itr = itr+1;
        if out_freq > 0 && mod(itr-1,out_freq)==0            
            disp(['itr : ' num2str(itr,'%3d') ' res ' num2str(max(norm(r)./normb),'%3.3e')]);
        end
        if itr >= maxit, break; end
        

    end
    
    
end
