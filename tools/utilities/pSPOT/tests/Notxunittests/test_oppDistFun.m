function test_suite = test_oppDistFun
% Experimental stage of oppDistFun
%test_oppCompositeFun  Unit tests for the oppCompositeFun operator
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function test_oppDistFun_fun
%% Testing oppDistFun with a fun
% Solving the problem of indexing over last dimension

m  = 500;
n  = 300;
o  = 5;
A1 = distributed.randn(m,n,o);
A2 = distributed.randn(m,o);
x  = distributed.randn(n,o);
F  = @pSPOT.test.funfun;
Q  = oppDistFun(A1,A2,F);
x  = x(:);
y = Q*x;

end %

function test_oppDistFun_numBlockDiag
%% Testing a oppDistFun version of oppNumBlockDiag
m = 500;
n = 300;
o = 5;
A = distributed.randn(m,n,o);
x = distributed.randn(n,o);
F = @pSPOT.test.funBlockDiag;
Q = oppDistFun(A,F);
B = oppNumBlockDiag(A);
x = x(:);

% Multiply
% fprintf('oppDistFun      :'); tic
y1 = Q*x; % toc
% fprintf('oppNumBlockDiag :'); tic
y2 = B*x; % toc

assertEqual(norm(y1-y2),0);

end



%% Ingenious code that will forever solve the problem of last-dimension
% indexing, from
% http://www.mathworks.com/matlabcentral/answers/846-multidimensional-colon-operator
% idx(1:ndims(A) - 1) = {':'};
% A(idx{:},t) = x;
%
