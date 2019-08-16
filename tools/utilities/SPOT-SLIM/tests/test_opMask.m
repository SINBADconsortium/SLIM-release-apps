function test_suite = test_opMask
%test_opMatrix  Unit tests for opMatrix.
initTestSuite;
end

function test_opMask_multiply
   %% Set up matrices and operators for problems
   M = opMask(true(3,1));
   x = randn(3) + 1i*randn(3);

   % Check opMask
   assertEqual(x,M*x);
   assertEqual(x',M*x');

end