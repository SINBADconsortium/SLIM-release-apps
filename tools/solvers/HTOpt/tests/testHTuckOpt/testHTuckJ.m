function test_suite = testHTuckJ
%TESTHTUCKJ - Ensures that the HT Jacobian + its adjoint behave as expected
%(O(h^2) Taylor error + correct adjoint test).
test_suite=buildFunctionHandleTestSuite(localfunctions);

    function testForwardError
        % Ensures that the Taylor error of the HT function decays as
        % O(h^2), averaged over 20 randomized trials.
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
        
        numTrials = 20;
        
        rateOfConvergence = zeros(numTrials,1);
        
        for j=1:numTrials
            x = project(dimTree.randn(),dimTree);
            J = opHTuckJ(dimTree,x);                        
            
            h = 10.^(-10:-3);
            X0 = dimTree.full(x);
            dx = project_horizontal(x,dimTree.randn(),dimTree);
            
            e = zeros(length(h),1);
            for i=1:length(e)
                e(i) = norm(dimTree.full(x + h(i)*dx) - X0 - h(i)*J*dx) ;
            end
                        
            rateOfConvergence(j) = log10(mean(e(2:end)./e(1:end-1)));
        end
        assertTrue(mean(abs(rateOfConvergence - 2)) < 0.2);
        
   function testAdjoint
       % Adjoint test for Jacobian
        n = 10; d = 5; dims = n*ones(1,d);
        kint = 10; kleaf = 5;
        dimTree = dimensionTree(dims,kleaf,kint);
                
        x = project(dimTree.randn(),dimTree);
        J = opHTuckJ(dimTree,x);
        
        dx = project_horizontal(x,randn(size(J,2),1),dimTree); dy = randn(size(J,1),1);
        s = J * dx; s = dy' * s;
        t = J' * dy; t = dx' * t;
        assertVectorsAlmostEqual(s,t,'relative');
        
        
        

