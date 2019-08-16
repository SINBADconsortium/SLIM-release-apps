function test_ToolsSolversLinearizedBregman
% Test for the lbm.m function

% Author: Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: May, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

fprintf('Starting linearized Bregman test run. \n');
tic

% dimensions of underdetermined system
m = 2000;
n = 4000;
A = randn(m,n);

% sparse solution solution with 1% non-zero entries
x = diag(rand(n,1)<=0.01)*randn(n,1); 

% measurements
bTrue = A*x;
% add noise
b = bTrue + 0.5*rand(size(bTrue));
sigma = norm(bTrue-b,2);


%% Set-Up Optimization
options.verbosity = 0;
options.iterations = 200;
options.xtrue = x;
x0 = [];
lambda = [];
 
% linearized Bregman with dynamic step size (default)
options.method = 's';
options.step = 'd';
[x1, res1] = lbm(A,b,lambda,sigma,x0,options);

% linearized Bregman with exact step size
options.method = 's';
options.step = 'e';
[x2, res2] = lbm(A,b,lambda,sigma,x0,options);

% Accelerated linearized Bregman with dynamic step size
options.method = 'a';
options.step = 'd';
[x3, res3] = lbm(A,b,lambda,sigma,x0,options);

% Accelerated linearized Bregman with exact step size
options.method = 'a';
options.step = 'e';
[x4, res4] = lbm(A,b,lambda,sigma,x0,options);

% error
err(1) = norm(x-x1,2);
err(2) = norm(x-x2,2);
err(3) = norm(x-x3,2);
err(4) = norm(x-x4,2);

% residual
res(1) = norm(res1,2);
res(2) = norm(res2,2);
res(3) = norm(res3,2);
res(4) = norm(res4,2);

for i=1:4
    % Error two-norm is smaller than 5 percent of the solution's two-norm
    assert(err(i) < 0.05*norm(x,2));
    % Residual is up to 10 percent within noise level
    assert(res(i) - sigma < 0.1*sigma);
end
t = toc;

fprintf(['PASSED in ',num2str(t),' seconds. \n']);

end
