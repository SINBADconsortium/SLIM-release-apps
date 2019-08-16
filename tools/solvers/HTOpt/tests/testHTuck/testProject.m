function testProject
%TESTPROJECT - Ensures that HT parameter orthogonalization doesn't change
%the tensor itself.
n = 10; d = 6; dims = n*ones(1,d);
kint = 10; kleaf = 5;
dimTree = dimensionTree(dims,kleaf,kint);

x = dimTree.randn();
x_orth = project(x,dimTree);

assertVectorsAlmostEqual(dimTree.full(x), dimTree.full(x_orth),'relative');

end

