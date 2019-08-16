function [U,S,V] = randsvd(A, rank,tol,maxrank)
% [U,S,V] = randsvd(A,rank,tol)
%    Randomized SVD for low-rank matrices based on Joel Tropp's
%    'Finding Structure with Randomness'
%
% Inputs
%    A    - low-rank matrix
%    rank - initial rank estimate
%    tol  - subspace error tolerance (default: 1e-2)
    if nargin < 2
        error('Need initial rank estimate');
    end
    if(exist('tol','var')==0)
        tol = 1e-2;
    end
	if(exist('maxrank','var')==0)
	        maxrank = size(A,1);
	end
    m = size(A,1);
    n = size(A,2);
    G = randn(n,rank + 5);
    R = A * G;
    [U,S,V] = svd(R,'econ');
    Q = U(:,1:min(rank,size(U,2)));
    err = norm(R - Q* (Q' * R),'fro');
    while(err > tol * S(1,1) & rank < maxrank)    
        rank = rank + round(0.005 * m);
        
        G = randn(n,round(0.005 * m));
        R = [R, A * G];
        [U,S,V] = svd(R,'econ');
        Q = U(:,1:min(rank,size(U,2)));
        err = norm(R - Q* (Q' * R),'fro');    
        disp(['Current rank ' num2str(rank) ' err ' num2str(err/S(1,1))]);
    end
    T = apply_op(A',Q);
    [U,S,V] = svd(T','econ');
    U = Q * U;

end

function W = apply_op(A,x)
    if(~isobject(A))
        W = A *x;
    else
        W = zeros(size(A,1),size(x,2));
        for i=1:size(x,2)
            W(:,i) = A * x(:,i);
        end
    end
        
end
