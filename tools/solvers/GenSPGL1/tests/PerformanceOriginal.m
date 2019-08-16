m = 1200; n = 5120; k = 200; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0 + 0.005 * randn(m,1);

sigma = 1e-3; 

[xOrig,rOrig,g,infoOrig] = spgl1Orig(A, b, 0, sigma); % Find BP sol'n.

[xNew,rNew,g,infoNew] = spgl1(A, b, 0, sigma); % Find BP sol'n.


opts = spgSetParms('funPenalty', @funLS2);

[xSq,rSq,g,infoSq] = spgl1(A, b, 0, sigma^2/2, [], opts); % Find BP sol'n.