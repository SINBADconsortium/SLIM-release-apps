function x = jacobi(A,alpha,b,x,mode,maxiter)
% Simply Jacobi overrelaxation smoother
%
% Curt Da Silva, 2016
%
% Usage:
%   y = jacobi(A,alpha,b,x,mode,maxiter);
%
% Input:
%   A     - matrix
%   b     - right hand side
%   x     - current estimate of the solution
%   alpha - overrelaxation parameter (default: 1.5)
%
% Output: 
%   y     - Jacobi-smoothed solution estimate
%
    if exist('alpha','var')==0||isempty(alpha), alpha = 1.5; end
    if mode~=1, A = A'; end
    d = diag(A);
    for i=1:maxiter
        x = x + alpha*(b-A*x)./d;
    end

end