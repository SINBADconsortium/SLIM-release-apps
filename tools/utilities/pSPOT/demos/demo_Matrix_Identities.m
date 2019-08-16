%   This page is for demoing useful matrix identities using pSpot
%   operators
%

% Kronecker
%%   kron(B',A)*vec(X) = vec(A*X*B) 
A = randn(5,10);
B = randn(20,10);
X = randn(10,20);
K = oppKron2Lo(B',A); 
C1 = A*X*B;
C1 = C1(:);
C2 = K*X(:); 
norm(C1 - C2) %% Answer should be 0

%% oppSweep(A) == oppKron2Lo(opDirac(ncols),A)
% Proof that oppSweep is equivalent to the kron of a dirac, but faster
n = 100;
A = opDFT(n);
x = distributed.randn(n,n,n);

S = oppSweep(A);
K = oppKron2Lo(opDirac(n),opKron(opDirac(n),A));

fprintf('C1 = S*x    : ');
tic, C1 = S*x; toc

fprintf('C2 = K*x(:) : ');
tic, C2 = K*x(:); toc % oppSweep is about 2x faster

norm( C1(:) - C2 ) % Answer should be 0

