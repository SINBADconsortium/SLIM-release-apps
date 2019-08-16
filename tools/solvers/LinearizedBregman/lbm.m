function [x, r, info] = lbm(A,b,lambda,sigma,x0,options)
% LBM Solve regularized basis pursuit or RBP denoise via the linearized
% Bregman method.
%
% minimize   lambda*||x||_1 + 0.5*||x||^2_2
% subject to ||A*x-b|| <= sigma
%
%--------------------------------------------------------------------------
%
%  [x, r, info] = lbm(A,b,lambda,sigma,x0,options)
%
% 
% INPUTS
% ======
% A        is an explicit n by m matrix or SPOT operator.
% b        is an n-size vector.
% lambda   is the parameter that weights the strongly convex regularizer in
%          the objective function. For lambda=[], lambda is determined from
%          the first application of A' to the residual and is set to the
%          absolute maximum value.
% sigma    is the noise level of the input data. For sigma=0, LB will try
%          to fit Ax=b exactly. For sigma != 0, LB will solve the basis
%          pursuit denoise problem (for lambda large enough).
% x0       is an m-size vector of the initial guess for the solution. For
%          x=[], LB will set x=0.
% options  is a structure for options. For options=[], LB will use all
%          default values.
%
% OUTPUTS
% =======
% x        is an m-size vector with the solution of the problem.
% r        is the residual r=Ax-b
% info     is a structure with the following information:
%          .lambda     default lambda value if no lambda was supplied by
%                      the user. 
%          .residual   Two-norm of the residual as a function of iterations
%          .error      Two-norm of the error, if the true solution was
%                      supplied by the user. Otherwise info.error=[]
%          .iterations final number of iterations.
%          .time       solution time [s]
%
% OPTIONS
% =======
% options.iterations   maximum number of iterations (default is 10).
%        .xtrue        true solution (if available) for calculating the 
%                      error two norm.
%        .step         'd' for dynamic step length (default), 'e' for exact
%                      step length using a line search
%        .method       's' for standard LB method (default), 'a' for
%                      accelerated LB
%        .verbosity    if '1', print information in each iteration
%

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


% Check options and set default values
options = checkOptions(options);

% default value for sigma
if isempty(sigma)
    sigma = 0;
end

% Initialize solution and variables
[m,n] = size(A);
if(isempty(x0))
    x = zeros(n,1);
else
    x = x0;
end
z = zeros(n,1);
v = zeros(n,1);
residual = ones(options.iterations,1)*sigma;
error = zeros(options.iterations,1);


% if solution not available return empty error array
if isempty(options.xtrue)
    error = [];
end

% initial residual and gradient
r = A*x-b;
rT = A'*r;

% heuristic value for lambda if not supplied
if(isempty(lambda))
    lambda = 1*max(abs((norm(r,2)/norm(rT,2))^2 * ...
             rT*max(0,1 - sigma/norm(r,2))));
end

% Main loop
iter = 0;
t = 1; alpha = 1; tic;
if options.verbosity
    fprintf('Iteration ||Residual||  ||Gradient||    Steplength\n'); 
end
while (norm(r,2)>sigma && iter < options.iterations)
    
    % calculate step length
    switch options.step
        case 'd'
            % dynamic step length
            t = (norm(r,2)/norm(rT,2))^2;
        case 'e'
            % exact step length
            tIni = t;
            fh = @(t) f1dim(t,z,rT,rT'*x  - norm(r,2)^2,lambda);
            t = fminsearch(fh,tIni);
    end
    
    switch options.method
        case 's'
            % gradient descent on dual variable
            z = z-t*rT*max(0,1 - sigma/norm(r,2));
            % update primal variable
            x = softThresholding(z,lambda);
        
        case 'a'
            % accelerated gradient descent
            vnew = z-t*rT*max(0,1 - sigma/norm(r,2));
            z = alpha*vnew+(1-alpha)*v;
            v = vnew;
            alpha = 1+iter/(iter+3);
            % update primal variable
            x = softThresholding(z,lambda);
    end
    
    % update residual and gradient
    r = A*x-b;
    rT = A'*r;
            
    iter = iter+1;
    if options.verbosity
        fprintf('    %d       %3.2e      %3.2e      %3.2e  \n',...
                iter, norm(r,2), norm(rT,2), t);
    end
    residual(iter) = norm(r,2);
    if (~isempty(options.xtrue))
        error(iter) = norm(x-options.xtrue,2);
    end
end

t=toc;
info.lambda = lambda;
info.sigma = sigma;
info.error = error;
info.residual = residual;
info.time = t;
info.iterations = iter;
if options.verbosity
    fprintf('Products with A : %d \n', iter+1);
    fprintf('Products with A'': %d \n', iter+1);
    fprintf('Total time [s]  : %2.2f \n\n', t);
end

end


function xSparse = softThresholding(x,lambda)
% Soft thresholding function
%
xSparse = max(abs(x)-lambda,0).*sign(x);
    
end


function res= f1dim(t, x, alpha, beta, lambda)
% F1DIM One-dimensional objective function for optimal step length
%
% res= f1dim(t, x, alpha, beta, lambda)
%
res = 0.5*norm(softThresholding(x-t*alpha,lambda),2)^2 + t*beta;
    
end


function options = checkOptions(options)
% CHECKOPTIONS Set default options if not supplied
%

% 10 iterations as default
if ~isfield(options,'iterations')
    options.iterations = 10;
end

% gradient descent as default
if ~isfield(options,'method')
    options.method = 's';
end

% dynamic step length as default
if ~isfield(options,'step')
    options.step = 'd';
end

% empty array as default
if ~isfield(options,'xtrue')
    options.xtrue = [];
end

% empty array as default
if ~isfield(options,'verbosity')
    options.verbosity = 1;
end

end


