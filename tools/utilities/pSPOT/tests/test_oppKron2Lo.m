function test_suite = test_oppKron2Lo
    initTestSuite;
end
    
function test_oppKron2Lo_builtin
%% Built-in unit tests for oppKron2Lo
m = 3;
n = 3;
A = opDCT(m);
B = opDFT(n);
K = oppKron2Lo(B,A);
utest(K,1);
end % builtin

function test_oppKron2Lo_basic
%% test_opKron  Unit tests for Kronecker products
   m = 2;
   n = 3;
   o = 4;
   A1 = randn(n,o) + 1i*randn(n,o);
   A2 = randn(n,m) + 1i*randn(n,m);
   A3 = randn(m,m) + 1i*randn(m,m);
   A  = kron(A1,kron(A2,A3));
   BB = opKron(opMatrix(A2),opMatrix(A3));
   B  = oppKron2Lo(opMatrix(A1),BB,1);
   x  = randn(size(A,1),1) + 1i*randn(size(A,1),1);
   y  = randn(size(A,2),1) + 1i*randn(size(A,2),1);
   assertElementsAlmostEqual(A *y, B *y)
   assertElementsAlmostEqual(A'*x, B'*x)
   assertElementsAlmostEqual(A ,double(B,1))
   Bleh = double(B,1);
   assertElementsAlmostEqual(A',Bleh')
end

function test_oppKron2Lo_emptylabs
%% Test for empty labs
% Setup x
spmd
    x = codistributed.randn(100,1,1);
    xpart = [1 zeros(1,numlabs-1)];
    xgsize = [100 1 1];
    xcodist = codistributor1d(3,xpart,xgsize);
    x = redistribute(x,xcodist);
end

A = opDFT(100);
K = oppKron2Lo(opDirac(1),A);
xvec = x(:);
y = K*xvec;
end % empty labs

function test_oppKron2Lo_dirac
%% Test for Dirac-skipping functionality
m = 3;
n = 3;
A = opDirac(m);
B = opDFT(n);
x = randn(n,m);
spmd
    x = codistributed(x,codistributor1d(2));
end
x = x(:);
A2 = A;
A2.isDirac = false;

% Dirac-Skipping
K1 = oppKron2Lo(A,B,1);

% Non-Skipping
K2 = oppKron2Lo(A2,B,1);

% fprintf('Dirac-skipping: '); tic; 
y1 = K1*x; 
% toc;

% fprintf('Non-skipping  : '); tic; 
y2 = K2*x; 
% toc;

assertEqual(y1,y2);
end % dirac

function test_oppKron2Lo_dirac_special
%% Dirac special
% Strange case encountered by Tristan
A  = randn(10,51);
K1 = opKron(opDirac(4),opKron(A,opDirac(101)));
K2 = oppKron2Lo(opDirac(4),opKron(A,opDirac(101)),1);
x1 = randn(5151,4);
x2 = distributed(x1);
y1 = K1*x1(:);
y2 = K2*x2(:);

assertElementsAlmostEqual(y1,y2);

end % dirac special

function test_oppKron2Lo_FoG
%% FoG
m  = 3;
A  = opDFT(m);
B  = opDFT(m*m);
K1 = B*opKron(A,A)*B;
K2 = B*oppKron2Lo(A,A,1)*B;
x  = K1.drandn;
assertElementsAlmostEqual(K1*x, K2*x);
end % FoG