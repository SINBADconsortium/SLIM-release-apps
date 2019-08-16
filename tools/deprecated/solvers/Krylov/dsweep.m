function x = dsweep(A,w,x,b,nsweeps)
% DSWEEP - CARP sweeps for preconditioning the linear system Ax = b
%
% Usage:
%    x = dsweep(A,w,x,b,nsweeps,n_threads);
% 
% Input:
%    A         - SPOT operator that supports carp sweeps OR a sparse matrix
%    w         - overrelaxation parameter 
%    x         - current iterate 
%    b         - right hand side
%    nsweeps   - number of sweeps
%    n_threads - number of threads
% 
% Output:
%    x         - CARP sweeps applied to x
%
    x = zeros(size(A,2),size(b,2));
    for i=1:size(A,2)
        c = w*(b(i,:) - A(i,:)*x)./norm(A(i,:))^2;
        x = x + (A(i,:)') * c;
    end
    for i=size(A,2):-1:1
        c = w*(b(i,:) - A(i,:)*x)./norm(A(i,:))^2;
        x = x + A(i,:)' * c;
    end
end