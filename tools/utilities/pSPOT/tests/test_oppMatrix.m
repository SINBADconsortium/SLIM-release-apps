function test_suite = test_oppMatrix
% test_oppMatrix  Unit tests for the oppMatrix operator
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_oppMatrix_builtin
%% Built-in tests for oppMatrix
m = 3; n = 3;
a = distributed.randn(m,n);
A = oppMatrix(a);
utest(A,3);

end % builtin

function test_oppMatrix_multiply
%% Multiplication test for oppMatrix
% There is a minor error between undistributed and distributed matrix
% multiplications because of the distributed nature of the distributed
% matrices. Error is in the order of e-15 so it is ok. We will just use
% assertElementsAlmostEqual
m  = 3; n = 3;
A1 = randn(m,n);
A2 = oppMatrix(distributed(A1));
x1 = [A2.drandn A2.drandn];
y1 = A1*x1;
y2 = A2*x1;
if isdistributed(y1), y1 = gather(y1); end
if isdistributed(y2), y2 = gather(y2); end
assertElementsAlmostEqual(y1,y2);

x2 = [A2.rrandn A2.rrandn];
z1 = A1'*x2;
z2 = A2'*x2;
if isdistributed(z1), z1 = gather(z1); end
if isdistributed(z2), z2 = gather(z2); end

assertElementsAlmostEqual(z1,z2);

end % multiply

function test_oppMatrix_basis_rn
%% Test for oppMatrix so that it generates the correct result
m  = 3; n = 3;
A1 = randn(m,n);
A2 = oppMatrix(distributed(A1));
x  = eye(n);
z  = A2*x;
if isdistributed(A1), A1 = gather(A1); end
if isdistributed(z),  z  = gather(z); end
assertEqual(z,A1);
end % basis

function test_oppMatrix_divide
%% test for divide of oppMatrix
% distributed matrix only supports left divide of square matrices
n  = 3;
A1 = randn(n);
A2 = oppMatrix(distributed(A1));
x  = A2.drandn;
y  = A1*x;
z  = A2\y;
if isdistributed(x), x = gather(x); end
if isdistributed(z), z = gather(z); end
assertElementsAlmostEqual(z,x);


end % divide

function test_oppMatrix_plus
%% test for plus of oppMatrix
m  = 3; n = 3;
A1 = randn(m,n);
A2 = oppMatrix(distributed(A1));
B  = randn(m,n);
y1 = A1 + B;
y2 = A2 + B;
z  = zeros(m,n);
assertEqual(z,gather(y1-double(y2)));
end % plus

function test_oppMatrix_minus
%% test for minus of oppMatrix
m  = 3; n = 3;
A1 = randn(m,n);
A2 = oppMatrix(distributed(A1));
B  = randn(m,n);
y1 = A1 - B;
y2 = A2 - B;
z  = zeros(m,n);
q  = y1 - double(y2);
if isdistributed(q), q = gather(q); end
assertEqual(z,q);
end % minus