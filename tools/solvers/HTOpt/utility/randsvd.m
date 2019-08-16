function [U,S,V] = randsvd(A, rank,tol)
% RANDSVD - Randomized SVD for computing SVDs of large scale,
% low-rank matrices. See 'Finding structure with randomness:
% probabilistic algorithms for constructing approximate matrix
% decompositions' - N. Halko, P.G. Matrinsson, J. A. Tropp for more details.
%   
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%  [U,S,V] = randsvd(A,rank,{tol});
%
% Input:
%  A       - large, low-rank matrix (if max(size(A)) < 1e3, uses
%            deterministic SVD)
%  rank    - estimated rank
%  tol     - relative error tolerance (default: 1e-2)
% 
% Output: 
%  U,S,V   - U - left singular vectors, V right singular vectors, 
%            S - diagonal matrix of singular values (same as svd)
%
if(exist('tol','var')==0)
    tol = 1e-2;
end
m = size(A,1);
n = size(A,2);
if m < 1e3 && n < 1e3
    [U,S,V] = svd(A,'econ');
else

G = randn(n,rank + 5);
R = A * G;
[U,S,V] = svd(R,'econ');
Q = U(:,1:min(rank,size(U,2)));
err = norm(R - Q* Q' * R,'fro');
while(err > tol * S(1,1))    
    rank = rank + round(0.1 * m);

    G = randn(n,round(0.1 * m));
    R = [R, A * G];
    [U,S,V] = svd(R,'econ');
    Q = U(:,1:min(rank,size(U,2)));
    err = norm(R - Q* Q' * R,'fro');    
    disp(['Current rank ' num2str(rank) ' rel err ' num2str(err/S(1,1))]);
end
T = apply_op(A',Q);
[U,S,V] = svd(T','econ');
U = Q * U;

end

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