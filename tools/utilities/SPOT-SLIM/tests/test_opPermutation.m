function test_suite = test_opPermutation
%test_opPermutation  Unit tests for opPermutation.
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function test_opPermutation_multiply

   rng('default');

   for i=1:2
   % Set up matrices and operators for problems
   p  = randperm(5)';
   A  = opPermutation(p);
   xr = randn(5,i);
   xi = sqrt(-1) * randn(5,i);
   x  = xr + xi;

   % Check opPermutation
   assertEqual(A * x  , x(p,:));
   assertEqual(A * xr , xr(p,:));
   assertEqual(A * xi , xi(p,:));

   assertEqual(A' * x(p,:)  , x);
   assertEqual(A' * xr(p,:) , xr);
   assertEqual(A' * xi(p,:) , xi);
   end

end

function test_opPermutation_divide

   rng('default');

   % Set up matrices and operators for problems
   for i=1:2
   p  = randperm(5)';
   A  = opPermutation(p);
   xr = randn(5,i);
   xi = sqrt(-1) * randn(5,i);
   x  = xr + xi;

   % Check opPermutation
   assertEqual(A \ x(p,:)  , x);
   assertEqual(A \ xr(p,:) , xr);
   assertEqual(A \ xi(p,:) , xi);

   assertEqual(A' \ x  , x(p,:));
   assertEqual(A' \ xr , xr(p,:));
   assertEqual(A' \ xi , xi(p,:));
   end

end
