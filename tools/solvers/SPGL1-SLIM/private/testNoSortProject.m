function test_suite = testNoSortProject
test_suite=buildFunctionHandleTestSuite(localfunctions);

    % BEGIN INITIALIZATION

    function opt = setup
        % Make sure environment is pristine
        clear all
        clear functions

        
        % Fix random generator for repeatable experiments
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',8888));
        
        opt.savedState = savedState;
        
    function teardown(opt)
        % Restore random stream and restore variables
        defaultStream.State = opt.savedState;
        clear opt
        
        % Remove the test SPGL1 from path
        clear all
        clear functions
    


    % BEGIN TESTS
    
    function testCorrectOutput(opt)
        for k = 1:5
            x = randn(1000000,1);
            x = abs(x);
            tau = rand(1)*norm(x,1);
            
            tic; x1 = oneProjectorMex(x,tau); toc;
            tic; x2 = oneProjectorMex_noSort(x,tau); toc;
            % norm(x1-x2,1)
            assertElementsAlmostEqual(x1, x2);
        end

    function testZeroTau(opt)
        x = randn(1000000,1);
        x = abs(x);
        tau = 0;
    
        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);
    
        assertElementsAlmostEqual(x1, x2);
        
    function testAlmostZeroTau(opt)
        x = randn(1000000,1);
        x = abs(x);
        tau = 0.0005;

        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);

        assertElementsAlmostEqual(x1, x2);

    function testLargeTau(opt)
        x = randn(1000000,1);
        x = abs(x);
        tau = 100*norm(x,1);

        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);

        assertElementsAlmostEqual(x1, x2);
        
    function testZeroVector(opt)
        x = zeros(1000000,1);
        x = abs(x);
        
        tau = 100;
        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);
        assertElementsAlmostEqual(x1, x2);
        
        tau = 0;
        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);
        assertElementsAlmostEqual(x1, x2);
        
    function testFeasibleVector(opt)
        x = randn(1000000,1);
        x = abs(x);
        tau = norm(x,1)+10;
        x(x < 0.5) = 0;
    
        x1 = oneProjectorMex(x,tau);
        x2 = oneProjectorMex_noSort(x,tau);
    
        assertElementsAlmostEqual(x1, x2);
        assertElementsAlmostEqual(x, x2);
        
