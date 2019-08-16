function testMatricize
%TESTMATRICIZE - ensures that matricization places tensor entries in the proper
%locations of the matrix

%Size of problem
n = 10; d = 4; 
dims = n*ones(1,d);
X = randn(dims);

%Singleton matricizations
for l=1:length(dims)
    compl_i = [1:l-1 l+1:d];
    compl_dims = dims(compl_i);
    X_i = matricize(X,l,compl_i);
    for i=1:dims(l)
        for j=1:prod(compl_dims)
            [j1,j2,j3] = ind2sub(compl_dims,j);
            J = zeros(d,1);
            J(compl_i) = [j1,j2,j3];
            J(l) = i;
            assertVectorsAlmostEqual(X(J(1),J(2),J(3),J(4)),X_i(i,j));
        end
    end
end

%Multiple dimensions
rows = [1 2]; cols = [3 4];
Xmat = matricize(X,rows,cols);
for i=1:prod(dims(rows))
   for j=1:prod(dims(cols))
       [i1,i2] = ind2sub(dims(rows),i);
       [j1,j2] = ind2sub(dims(cols),j);
       assertVectorsAlmostEqual(X(i1,i2,j1,j2),Xmat(i,j));
   end
end

%Permutation of dimensions
rows = [2 1]; cols = [3 4];
Xmat = matricize(X,rows,cols);
for i=1:prod(dims(rows))
   for j=1:prod(dims(cols))
       [i1,i2] = ind2sub(dims(rows),i);
       [j1,j2] = ind2sub(dims(cols),j);
       assertVectorsAlmostEqual(X(i2,i1,j1,j2),Xmat(i,j));
   end
end

end

