function testDematricize
%TESTMATRICIZE - ensures that matricization places tensor entries in the proper
%locations of the matrix
n = 10; d = 4; 
dims = n*ones(1,d);
X = randn(dims);

%Singleton matricizations
for i=1:length(dims)
    compl_i = [1:i-1 i+1:d];
    Xmat = matricize(X,i,compl_i);
    Xdemat = dematricize(Xmat,dims,i,compl_i);
    assertVectorsAlmostEqual(vec(X),vec(Xdemat),'relative');
end

rows = [1 2]; cols = [3 4];
Xmat = matricize(X,rows,cols);
Xdemat = dematricize(Xmat,dims,rows,cols);
assertVectorsAlmostEqual(vec(X),vec(Xdemat),'relative');

rows = [1 2]; cols = [4 3];
Xmat = matricize(X,rows,cols);
Xdemat = dematricize(Xmat,dims,rows,cols);
assertVectorsAlmostEqual(vec(X),vec(Xdemat),'relative');

rows = [2 1]; cols = [3 4];
Xmat = matricize(X,rows,cols);
Xdemat = dematricize(Xmat,dims,rows,cols);
assertVectorsAlmostEqual(vec(X),vec(Xdemat),'relative');
end

