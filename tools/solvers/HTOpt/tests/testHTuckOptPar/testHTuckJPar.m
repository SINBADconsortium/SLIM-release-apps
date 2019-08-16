function testHTuckJPar
% This function ensures that opHTuckJ and oppHTuckJ compute the
% same thing.

n = 10; d = 5; dims = n*ones(1,d);
kint = 10; kleaf = 5;
dimTree = dimensionTree(dims,kleaf,kint);

x = project(dimTree.randn(),dimTree);

X = dimTree.full(x);
Xdist = dimTree.fullDist(x);

spmd
    codist = getCodistributor(Xdist);
    Xdistloc = getLocalPart(Xdist);
    rloc = Xdistloc;
    bloc = randn(size(rloc));
    rloc = rloc - bloc;
    r = codistributed.build(rloc,codist,'noCommunication');
end

rloc = gather(r);
J = opHTuckJ(dimTree,x);
Jpar = oppHTuckJ(dimTree,x);
dx = J' * rloc; dx2 = Jpar' * r;
assertVectorsAlmostEqual(dx,dx2);

J = opHTuckJ2(dimTree,x);
Jpar = oppHTuckJ2(dimTree,x);
dx = J' * rloc; dx2 = Jpar' * r;
assertVectorsAlmostEqual(dx,dx2);

end
