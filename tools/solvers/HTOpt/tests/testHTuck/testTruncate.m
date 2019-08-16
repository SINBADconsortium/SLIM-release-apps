function test_suite = testTruncate
    test_suite=buildFunctionHandleTestSuite(localfunctions);
    function testTruncate1
        %TESTTRUNCATE1 - ensures that truncating a HT tensor results in the same
        %tensor afterwards
        
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
        
        x = project(dimTree.randn(),dimTree);
        X = dimTree.fullND(x);
        y = dimTree.truncate(X);
        
        assertVectorsAlmostEqual(vec(X),dimTree.full(y),'relative');

    
    function testTruncate2
        %TESTTRUNCATE2 - ensures that truncating a general HT tensor
        % to a particular accuracy yields the same accuracy
        
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
        
        x = project(dimTree.randn(),dimTree);
        X = dimTree.fullND(x);
        
        [dimTree,y] = truncate_ht(X, 1e-6, true);
        assertVectorsAlmostEqual(vec(X),dimTree.full(y),'relative');
                
 

