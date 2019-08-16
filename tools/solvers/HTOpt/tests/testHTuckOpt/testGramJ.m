function test_suite = testGramJ
%TESTHTUCKJ - Ensures that the HT Gramian Jacobian + its adjoint behave as expected
%(O(h^2) Taylor error + correct adjoint test).
initTestSuite;

    function testForwardError
        % Ensures that the Taylor error of the HT function decays as
        % O(h^2), averaged over 20 randomized trials.
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
        
        numTrials = 10;
        
        rateOfConvergence = zeros(numTrials,1);
        
        for j=1:numTrials
            x = project(dimTree.randn(),dimTree);
            J = opGramianJ(dimTree,x);                        
            
            h = 10.^(-10:-3);
            X0 = dimTree.toVecGram(dimTree.gramian(x));
            dx = dimTree.randn();
            
            
            e = zeros(length(h),1);
            for i=1:length(e)
                e(i) = norm(dimTree.toVecGram(dimTree.gramian(x + h(i)*dx)) - X0 - h(i)*J*dx) ;
            end
                        
            rateOfConvergence(j) = log10(mean(e(2:end)./e(1:end-1)));
        end
        assertTrue(mean(abs(rateOfConvergence - 2)) < 0.2);
        
   function testAdjoint
       % Adjoint test for Gramian Jacobian
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
                
        x = project(dimTree.randn(),dimTree);
        J = opGramianJ(dimTree,x);
        
        dx = randn(size(J,2),1); dy = randn(size(J,1),1);
        s = J * dx; s = dy' * s;
        t = J' * dy; t = dx' * t;
        assertVectorsAlmostEqual(s,t,'relative');
        
        
        

