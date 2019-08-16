function testSparseObj
% Test that the dense + sparse objectives are giving the same result
    
    n = 10; d = 5; dims = n*ones(1,d);
    kint = 10; kleaf = 5;
    dimTree = dimensionTree(dims,kleaf,kint);
    
    x = project(dimTree.randn(),dimTree);
    e = false(prod(dims),1);
    e(randperm(prod(dims),100*n)) = 1;
    I = find(e);
    Ifull = idx1d2nd(dims,I);
    b = randn(prod(dims),1);
    X = dimTree.full(x);
    X(~e) = 0;
    bsub = b(I);
    b(~e) = 0;
    r = (X - b);
    f = 0.5*norm(r)^2;
    f2 = LSMisfitHT(x,Ifull,bsub,dimTree);
            
    % Objectives should be the same
    assertVectorsAlmostEqual(f,f2,'relative');
    
    J = opHTuckJ2(dimTree,x);
    g = J' * r;
    
    [f2,g2] = LSMisfitHT(x,Ifull,bsub,dimTree);
    
    % Objectives should be the same
    assertVectorsAlmostEqual(f,f2,'relative');
    
    % Objectives should be the same
    assertVectorsAlmostEqual(g,g2,'relative');

    
    
    
    
    
