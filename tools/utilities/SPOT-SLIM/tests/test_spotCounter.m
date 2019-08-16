function test_suite = test_spotCounter
%test_indexing  Unit tests for indexing
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function test_spotCounter_forward_adjoint
%%
A = opDirac(10);
x = ones(10,1);

% apply forward mode 3 times and adjoint mode 2 times
for k=1:3
    A*x;
end
for k=1:2
    A'*x;
end

% A.nprods should be [3 2]
assertEqual(A.nprods(1), 3)
assertEqual(A.nprods(2), 2)

end
