function [ok] = spgl1Test()

m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0 + 0.005 * randn(m,1);

% outliers
s = 5;
Err = zeros(m, 1);
pErr = randperm(m);
Err(pErr(1:s)) = randn(s,1);
b = b + Err;

opts = spgSetParms('optTol',1e-4, ...
                   'project', @NormL1_project, ...
                   'primal_norm', @NormL1_primal, ...
                   'dual_norm', @NormL1_dual, ...
                   'subspaceMin', 0 ...
                   );


params.nu = 1e-2;
params.hub = 0.005;
params.rho = 1;

%% Run three separate root finding problems

%sigma = 1e-1;
sigma = funLS(Err, params);


%params.funForward = @funForward;
params.funPenalty = @funLS;


[xLS, r, g, info] = spgl1(A, b, [], sigma, [], opts, params);
%[xLS, r, g, info] = spgl1Classic(A, b, [], sigma, 0*x0, opts);

%%

opts.funPenalty = @funHUB;
%sigma = 5e-2;
sigma = funHUB(Err, params);
[xHUB, r, g, info] = spgl1(A, b, [], sigma, [], opts, params );

opts.funPenalty = @funST;
%sigma = 1e-1;
sigma = funST(Err, params);
[xST, r, g, info] = spgl1(A, b, [], sigma, [], opts, params );


%sigma = 1e0;
%[xMY, r, g, info] = nlspgl1(@funForward, @funMY, params, b, [], sigma, 0*x0, opts );


%% Plot the results

figure(1)
plot(1:n, x0 + 2, 1:n, xLS(1:n), 1:n, xHUB(1:n) - 2, 1:n, xST(1:n) -4,  '-k', 'Linewidth', 1.5);

legend('true', 'LS', 'Hub', 'ST');

figure(2);
plot(1:m, b - A*x0 + 1.5, 1:m, b - A*xLS(1:n), 1:m, b - A*xHUB(1:n) - 1.5,  1:m, b-A*xST(1:n) - 3, '-k', 'Linewidth', 2); 
legend('true', 'LS', 'Hub', 'ST');


end

