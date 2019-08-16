function testSpan
%TESTSPAN - ensures that for the HT tensor X, 
%   span(X^{(t)}) = span(U_t) for each t \in T \setminus t_root
%
n = 10; d = 5; dims = n*ones(1,d);
kint = 10; kleaf = 5;
dimTree = dimensionTree(dims,kleaf,kint);

x = project(dimTree.randn(),dimTree);
X = dimTree.fullND(x);
fullTree = dimTree.fullTree(x);
itr = dimTree.iterator('down');

%Subspace projector
projector = @(A) A *  A';

%Subspace distance metric
dist = @(A,B) norm(projector(A) - projector(B));

while itr.advance()
   if ~itr.isRoot()
       [U_full,~,~] = svd(matricize(X,itr.getDims()));
       rank = itr.getRank();
       U_full = U_full(:,1:rank);
       U_ht = fullTree{itr.depth}{itr.idx}.U;
       assertTrue(dist(U_full,U_ht) < 1e-6);
   end
end

end