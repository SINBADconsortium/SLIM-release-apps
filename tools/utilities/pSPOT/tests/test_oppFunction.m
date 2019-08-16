function test_suite = test_oppFunction
%test_oppDictionary  Unit tests for the Dictionary meta operator
initTestSuite;
end

function test_oppFunction_zeros_preallocation
%% Test for zeros allocation class bug
x = distributed.randn(5,2);
F = oppFunction(5,5,@(x,mode) x);
assertEqual(x,F*x);
assertEqual(x,F'*x);
end