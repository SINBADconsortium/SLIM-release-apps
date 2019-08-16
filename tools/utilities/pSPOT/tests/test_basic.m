function test_suite = test_basic
%test_oppBlockDiag  Unit tests for the opBlockDiag operator
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function test_basic_mtimesdistribution
%% Test to assure the consistency of distribution after mtimes
    spmd
        x = codistributed.randn(10,1);
        part = zeros(1,numlabs);
        part(1) = 10;
        x = redistribute(x,codistributor1d(1,part));
    end
    
    y1 = 2*x;
    y1 = oppDirac(10)*y1;
    y1 = 0.1*y1;
    
    y2 = 0.5*oppDirac(10)*2*x;
    
    spmd
        l1 = length(getLocalPart(y1));
        l2 = length(getLocalPart(y1));
    end
    
    assertEqual([l1{:}],[l2{:}]);
end