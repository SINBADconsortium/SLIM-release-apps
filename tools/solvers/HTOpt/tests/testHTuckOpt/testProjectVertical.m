function testProjectVertical
%%TESTPROJECTVERTICAL - Ensures that project_vertical extracts the
%%component of null(opHTuckJ)
    n = 10; d = 5; dims = n*ones(1,d);
    kint = 10; kleaf = 5;
    dimTree = dimensionTree(dims,kleaf,kint);
    
    x = project(dimTree.randn(),dimTree);
    dx = dimTree.randn();  
    J = opHTuckJ(dimTree,x);
     
    %Filtering out nullspace    
    z = J * (project_vertical(x,dx,dimTree));
    assertVectorsAlmostEqual(norm(z),0,'absolute',1e-6);
    
    y = J * dx;
    z = J * (dx - project_vertical(x,dx,dimTree));
    assertVectorsAlmostEqual(y,z,'relative',1e-6); 
   
    
    
