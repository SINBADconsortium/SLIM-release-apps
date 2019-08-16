function testGaussNewton

    n = 10; d = 5; dims = n*ones(1,d);
    kint = 10; kleaf = 5;
    dimTree = dimensionTree(dims,kleaf,kint);
    
    x = project(dimTree.randn(),dimTree);
    dx = project_horizontal(x,dimTree.randn(),dimTree);    
    J = opHTuckJ(dimTree,x);
    H = opHTuckGN(dimTree,x);
    y = J' * J * dx; z = H * dx;
    
    %Hessian equivalence
    assertVectorsAlmostEqual(y,z,'relative');
    
    %Inverse Hessian
    assertVectorsAlmostEqual(dx,H\z,'relative');
    
    
    
