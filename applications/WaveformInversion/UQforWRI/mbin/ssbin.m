function x = ssbin(A,nmv,n)
% Stochastic matrix-free binormalization for symmetric real A.
% x = ssbin(A,nmv,n)
%   A is a symmetric real matrix or function handle. If it is a
%     function handle, then v = A(x) returns A*x.
%   nmv is the number of matrix-vector products to perform.
%   [n] is the size of the matrix. It is necessary to specify n
%     only if A is a function handle.
%   diag(x) A diag(x) is approximately binormalized.
  
% Jan 2010. Algorithm and code by Andrew M. Bradley (ambrad@stanford.edu).
% Aug 2010. Modified to record and use dp. New omega schedule after running
%   some tests.
% Jul 2011. New strategy to deal with reducible matrices: Use the old
%   iteration in the early iterations; then switch to snbin-like behavior,
%   which deals properly with oscillation.
  
  op = isa(A,'function_handle');
  if(~op) n = size(A,1); end
  d = ones(n,1); dp = d;
  for(k = 1:nmv)
    % Approximate matrix-vector product
    u = randn(n,1);
    s = u./sqrt(dp);
    if(op) y = A(s); else y = A*s; end
    % omega^k
    alpha = (k - 1)/nmv;
    omega = (1 - alpha)*1/2 + alpha*1/nmv;
    % Iteration
    d = (1-omega)*d/sum(d) + omega*y.^2/sum(y.^2);
    if (k < min(32,floor(nmv/2)))
      % First focus on making d a decent approximation
      dp = d;
    else
      % This block makes ssbin behave like snbin except for omega
      tmp = dp;
      dp = d;
      d = tmp;
    end
  end
  x = 1./(d.*dp).^(1/4);
