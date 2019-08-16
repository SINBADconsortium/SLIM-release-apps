function [dm,r,g,info] = L1Solver(J,dU,C,options)
%Solve the following Basis Pursuit problem:
%
%   minimize ||dm||_1 subject to Jdm = dU,
%
%   where J is the Jacobian (an M x N matrix), dU is the data residual (an M-vector),
%   and dm is the model residual (an N-vector). SIGMA is a non negative scalar
%
%   C is the curvelet operator
%
%   dm = L1Solver(J,dU,C) solves the problem
%
%   dm = L1Solver(J,dU,C,options) specifies options that are set using SPGSETPARMS
%
%   [dm,R,G,INFO] = L1Solver(J,dU,C,options) additionally returns the residual 
%   R = B - A*X (which should be small), the objective gradient 
%   G = A'*R, and an INFO structure.  (See SPGL1 for a description of this
%   last output argument.)   
%
%   See also spgl1, spgSetParms, spg_bpdn, spg_lasso.


%   Derived from spg_bp.m , Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1

if ~exist('options','var'), options = []; end
if ~exist('dU','var') || isempty(dU)
    error('Second argument cannot be empty.');
end
if ~exist('J','var') || isempty(J)
    error('First argument cannot be empty.');
end

sigma = 0;
tau = 0;
x0  = [];

if isfield(options,'tau')
   tau = options.tau;
end

if isfield(options,'sigma')
   sigma = options.sigma;
end

if isfield(options,'x0')
   x0 = C*options.x0; 
end


[dm,r,g,info] = spgl1(J,dU,tau,sigma,x0,options);
