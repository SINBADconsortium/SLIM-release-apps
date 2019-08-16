function [x,log]=Compute_projection_ADMM(x,A,funProj,options,R)
%ADMM function to solve a class of projection problems of a vector (y) onto a
%constraint set (C), possibly defined in a transform domain. Theory assumes
%a closed and convex set C. Algorithm may still work when a nonconvex set
%is used.
%
% Solves: min_(x,z) (1/2)||x-y||^2_2 s.t. Ax in C
% formulated as
% min_(x,z) (1/2)||x-y||^2_2 + I_C(z) s.t. Ax=z
% where I_C is the indicator function for set C and A is a matrix defining
% the transform domain (discrete gradient, TV, a basis, etc)
%
% This version uses the scaled form ADMM. Most of the code (and notation) is modeled
% after:
% Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction method of multipliers." Foundations and Trends® in Machine Learning 3.1 (2011): 1-122.
%
% input:
%       x                   -   vector to be projected onto C
%       funProj             -   funProj(input) projects input on C (output as vector)
%       A                   -   Transform domain operator, explicit (sparse) matrix or SPOT operator
%       options.
%           options.maxit   -   maximum number of iterations
%           options.evol_rel_tol- tolerance on relative evolution between iterations: exit if norm(x-x_old)/norm(x_old) becomes too small
%           options.rho     -   initial vaue of penalty parameter
%           options.adjust_rho - if (=1) adjust rho heuristically
%           options.feas_tol-  feasibility tolerance for warning message, not for stopping condition          
%       R                   -   (optional) Cholesky factor of (...) in case factorization caching is used in combination with a fixed augmented-Lagrangian penalty parameter rho
%
% output:
%       x           -   result
%       log.        -   constains log information about various quantities per iteration

% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date:January 2016.

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)

%% Parse options
maxit    = 500;
evol_rel_tol = 2e-5;
rho = 1e2;
adjust_rho = 1;
mu = 10;
tau_increment = 2;
tau_decrement = 2;
feas_tol = 5e-2;
nonconvex=0;
if exist('options','var') && ~isempty(options)
    if isfield(options,'maxit'),         maxit        = options.maxit; end
    if isfield(options,'evol_rel_tol'),  evol_rel_tol = options.evol_rel_tol; end
    if isfield(options,'rho'),           rho          = options.rho; end
    if isfield(options,'adjust_rho'),    adjust_rho   = options.adjust_rho; end
    if isfield(options,'mu'),            mu           = options.mu; end
    if isfield(options,'tau_increment'), tau_increment= options.tau_increment; end
    if isfield(options,'tau_decrement'), tau_decrement= options.tau_decrement; end
    if isfield(options,'feas_tol'),      feas_tol     = options.feas_tol; end
    if isfield(options,'nonconvex'),     nonconvex    = options.nonconvex; end
end

%% Input checks
% ADMM can work with imaginary numbers and zeros, but in geophysical
% inverse problems the input and output should typically be real and
% positive.
if isreal(x)==0
    disp('warning:input for ADMM is not real')
end
if nnz(x)~=length(x)
    disp('warning: input of ADMM contains ZEROS')
end

%% define stopping criteria constants
%hardcoded for now, change this based on the problem
ABSTOL   = 5e-7;
RELTOL   = 5e-6;
[pp,nn]  = size(A);

%% initialize
y=x;
z=A*x;
u=zeros(length(z),1);

for i=1:maxit
    x_old=x;
    %% x-minimization
    %either take in a cholesky (QR not implemented yet) that needs to be formed once every frequency batch
    %or solve iteratively
    if adjust_rho==1
        [x,flag,relres,iter,resvec] = lsqr([speye(length(x));sqrt(rho)*A],[y;sqrt(rho)*(z-u)],1e-6,1e3,[],[],x); %initial guess is the current x-estimate (warmstart)
        if iter==0; x=x_old; end;
    else
        x = R\(R'\(y+rho*A'*(z-u)));
    end
    
    %we need this mat-vec product multiple times, so do it once and save it
    Ax=A*x;
    
    %% z-minimization
    z_old = z;
    
    % projection onto the constraint set of Ax+u
    z = funProj(Ax+u);
    
    %% u-update
    u = u + Ax - z;
    
    %% logging
    r_pri_norm  = norm(Ax-z);
    r_dual_norm =  norm(rho*A'*(z - z_old));
    
    log.r_pri_norm(i)   = r_pri_norm;
    log.obj(i)          = (0.5*norm(x-y)^2);
    log.dual_norm(i)    = r_dual_norm;
    log.rho(i)          = rho;
    
    
    %% Stopping conditions
    log.eps_pri(i) = sqrt(pp)*ABSTOL + RELTOL*max(norm(Ax), norm(-z));
    log.eps_dual(i)= sqrt(nn)*ABSTOL + RELTOL*norm(rho*A'*u);

    if (r_pri_norm < log.eps_pri(i) && r_dual_norm < log.eps_dual(i))
        disp(['satisfies stopping conditions, exiting ADMM (iteration ',num2str(i),' )'])
        break;
    end
      
    if norm(x_old-x)/norm(x)<evol_rel_tol && i>3
        disp(['relative evolution to small, exiting ADMM (iteration ',num2str(i),' )'])
       break 
    end
    
    %additionaly, if the constraint set is non-convex:
    if nonconvex
        
    end
    
    %% rho (penalty parameter) adjustment (heuristic)
    if adjust_rho == 1
        if r_pri_norm > mu*r_dual_norm
            rho = rho * tau_increment;
            u   = u./tau_increment;
        elseif r_dual_norm>mu*r_pri_norm
            rho = rho./tau_decrement;
            u   = u.*tau_decrement;
        else
            %rho = rho; %do nothing
        end
    end
    
end

%% Output checks

%feasibility check
nAx      = norm(Ax);
feas_err = norm(funProj(Ax)-Ax)/nAx;
if feas_err>feas_tol && nAx>1000*eps
    disp(['iteration: ',num2str(i),',  final feasibility error ADMM : ',num2str(feas_err)])
    disp(['iteration: ',num2str(i),',  final primal feasibility error ADMM : ',num2str(r_pri_norm/norm(x))])
end
if isreal(x)==0
    disp('warning: Result of ADMM is not real')
end
if nnz(x)~=length(x)
    disp('warning: Result of ADMM contains ZEROS')
end
