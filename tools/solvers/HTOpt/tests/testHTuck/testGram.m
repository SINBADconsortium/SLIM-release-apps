function testGram
%TESTGRAM - Ensures that the singular values of the Gramian matrices
%correspond to the singular values of the matricizations of X
n = 10; d = 5; dims = n*ones(1,d);
kint = 10; kleaf = 5;
dimTree = dimensionTree(dims,kleaf,kint);

x = project(dimTree.randn(),dimTree);
X = dimTree.fullND(x);
G = dimTree.gramian(x);
T = dimTree.emptyTree(); 
T = dimTree.copyCell(T,G,'G');
itr = dimTree.iterator('down',T);

while itr.advance()
   if ~itr.isRoot()
       S = svd(matricize(X,itr.getDims()));
       rank = itr.getRank();
       S = S(1:rank);
       S_gram = svd(itr.getValue('G'));
       assertVectorsAlmostEqual(sqrt(S_gram),S,'relative');
   end
end

end

