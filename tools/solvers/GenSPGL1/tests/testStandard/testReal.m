function test_suite = testReal
test_suite=buildFunctionHandleTestSuite(localfunctions);

    % BEGIN INITIALIZATION

    function opt = setup
        % Make sure environment is pristine
        clear all
        clear functions
        
        % Make sure SPGL1 exists in path
        addpath('../..')
        
        % Fix random generator for repeatable experiments
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',8888));
        
        % Make test cases
        m = 300;
        n = 1000;
        k = 20;
        
        A = randn(m,n);
        x = zeros(n,1);
        x(floor((n-1)*rand(k,1))+1) = randn(k,1); % inject k-sparse signal in random places
        b = A*x;
        
        % Set sparse recovery element-wise relative tolerance
        recov_tol = 4e-6;
        
        % Inject options into test cases
        opt.A = A;
        opt.b = b;
        opt.x = x;
        opt.recov_tol = recov_tol;
        opt.savedState = savedState;

    function teardown(opt)
        % Restore random stream and restore variables
        defaultStream.State = opt.savedState;
        clear opt
        
        % Remove the test SPGL1 from path
        rmpath('../..')
        clear all
        clear functions
    


    % BEGIN TESTS
    
    function testRealSPGL1(opt)
        % figure;plot(opt.x)
        spg_opts  = spgSetParms('verbosity' ,1);
        [x, r, g, info] = spgl1(opt.A, opt.b, [], 1e-5, [], spg_opts);
        % figure;plot(x)
        % figure;plot(x-opt.x)
        norm(x-opt.x,inf)
        assertElementsAlmostEqual(x, opt.x, 'absolute', opt.recov_tol);
        assertTrue(info.iter < 70);
