function [x,log]=Compute_projection_ADMM_A_unscaled(x,A,funProj,options,ini_guess)
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
maxit           = 500;
evol_rel_tol    = 1e-5;
feas_tol        = 1e-2;
obj_tol         = 1e-3;
rho             = 1e0;
adjust_rho      = 1;
adjust_rho_type = 'BB';
adjust_gamma    = 1;
gamma           = 1;
mu              = 10;
tau_increment   = 2;
tau_decrement   = 2;

if exist('options','var') && ~isempty(options)
    if isfield(options,'maxit'),         maxit        = options.maxit; end
    if isfield(options,'evol_rel_tol'),  evol_rel_tol = options.evol_rel_tol; end
    if isfield(options,'rho'),           rho          = options.rho; end
    if isfield(options,'gamma'),         gamma        = options.gamma; end
    if isfield(options,'adjust_rho'),    adjust_rho   = options.adjust_rho; end
    if isfield(options,'adjust_gamma'),  adjust_gamma = options.adjust_gamma; end
    if isfield(options,'adjust_rho_type'),adjust_rho_type   = options.adjust_rho_type; end
    if isfield(options,'mu'),            mu           = options.mu; end
    if isfield(options,'tau_increment'), tau_increment= options.tau_increment; end
    if isfield(options,'tau_decrement'), tau_decrement= options.tau_decrement; end
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

%% initialize
counter=1;
ind_ref = inf;
m=x;
if nargin>4
    x=ini_guess;
end
z=A*x;
l=zeros(length(z),1);
l_hat_0=l;
l_0=l;
x_0=x;
z_0=z;

I = speye(length(x));
AA=A'*A;
for i=1:maxit
    x_old=x;
    %% x-minimization
%    [x,flag,relres,iter,resvec,lsres] = lsqr([I;sqrt(rho)*A],[m;sqrt(rho)*z+l/sqrt(rho)],1e-7,1e3,[],[],x); %initial guess is the current x-estimate (warmstart)
%    log.lsqr_relres(i)=lsres(end);
     Q=I+rho*AA;%[I;sqrt(rho)*A]'*[I;sqrt(rho)*A];
     rhs = [I;sqrt(rho)*A]'*[m;sqrt(rho)*z+l/sqrt(rho)];
    [x,flag,relres,iter,resvec] = pcg(Q,rhs,1e-7,1e3,[],[],x);
    log.lsqr_relres(i)=relres(end);
    log.lsqr_it(i)=iter;
   
    if iter==0; x=x_old; end;
    
    %we need this mat-vec product multiple times, so do it once and save it
    Ax=A*x;
    
    %% relaxation step
    x_hat = gamma*Ax + (1-gamma)*z;
    
    %% z-minimization
    z_old = z;
    % projection onto the constraint set of Ax+u
    z = funProj(x_hat-l/rho);
    
    %% u-update
    l_old = l;
    l = l + rho*(-x_hat + z);
%     if i==1
%         l_hat_old =  l_old + rho* ( -Ax + z_old);
%     end


    
    %% logging
    r_pri  = norm(Ax-z);
    r_dual = norm(rho*A'*(z - z_old));
    
    log.r_pri(i)        = r_pri;
    log.r_dual(i)       = r_dual;
    log.obj(i)          = (0.5*norm(x-m)^2);
    log.rho(i)          = rho;
    log.gamma(i)        = gamma;
    log.evol_x(i)       = norm(x_old-x)/norm(x);
    
    %log feasibility
    if mod(i,25)==0 %log every 25 it
        log.set_feasibility(counter)=norm(funProj(Ax)-Ax)/norm(Ax);
        counter=counter+1;
    end
    
    %% Stopping conditions
    [stop,adjust_rho,ind_ref]=stop_admm(log,i,evol_rel_tol,feas_tol,obj_tol,adjust_rho,ind_ref);
    if stop==1
        output_check_ADMM(Ax,funProj,x,i);
        break
    end
    
    %%adjust penalty parameter
    if adjust_rho == 1 || adjust_gamma == 1 || strcmp(adjust_rho_type,'RB')==1
        [rho,gamma,l_hat_0,l_0,x_0,z_0]=adapt_rho(i,gamma,rho,adjust_rho_type,adjust_gamma,adjust_rho,r_pri,r_dual,mu,tau_increment,tau_decrement,z,z_old,x,x_old,Ax,A,l,l_old,l_hat_0,l_0,x_0,z_0);
    end
end

output_check_ADMM(Ax,funProj,x,i);
end

%% penalty parameter adaptation function
function [rho,gamma,l_hat_0,l_0,x_0,z_0]=adapt_rho(i,gamma,rho,adjust_rho_type,adjust_gamma,adjust_rho,r_pri,r_dual,mu,tau_increment,tau_decrement,z,z_old,x,x_old,Ax,A,l,l_old,l_hat_0,l_0,x_0,z_0)
% if i>50
%     adjust_rho_type = 'BB';
% end
switch adjust_rho_type
    case 'RB'        % residual balancing
        if r_pri > mu*r_dual
            rho = rho * tau_increment;
        elseif r_dual>mu*r_pri
            rho = rho./tau_decrement;
        else
            %rho = rho; %do nothing
        end
    case 'BB' % Barzilai-Borwein type for Douglash-Rachford splitting on the dual problem
        eps_correlation = 0.2; %hardcored and suggested value by the paper based on numerical evidence
        l_hat     = l_old + rho* ( -Ax + z_old);%%cont here
        if mod(i,3)==0 && i>1
            l_hat     = l_old + rho* ( -Ax + z_old);%%cont here
            d_l_hat = l_hat - l_hat_0; 
            d_H_hat = A*(x-x_0);
            alpha_hat_MG = (d_H_hat'*d_l_hat)/(d_H_hat'*d_H_hat);
            alpha_hat_SD = (d_l_hat'*d_l_hat)/(d_H_hat'*d_l_hat);
            if (2*alpha_hat_MG) > alpha_hat_SD
                alpha_hat = alpha_hat_MG;
            else
                alpha_hat = alpha_hat_SD - alpha_hat_MG/2;
            end
            alpha_correlation =  (d_H_hat'*d_l_hat)/(norm(d_H_hat)*norm(d_l_hat));
            
            d_l = l - l_0;
            d_G_hat = -(z-z_0);
            beta_hat_MG = (d_G_hat'*d_l)/(d_G_hat'*d_G_hat);
            beta_hat_SD = (d_l'*d_l)/(d_G_hat'*d_l);
            if (2*beta_hat_MG) > beta_hat_SD
                beta_hat = beta_hat_MG;
            else
                beta_hat = beta_hat_SD - beta_hat_MG/2;
            end
            beta_correlation =  (d_G_hat'*d_l)/(norm(d_G_hat)*norm(d_l));
            
            %update rho
            if adjust_rho == 1
            if alpha_correlation > eps_correlation && beta_correlation > eps_correlation
                rho = sqrt(alpha_hat*beta_hat);
            elseif alpha_correlation > eps_correlation && beta_correlation <= eps_correlation
                rho = alpha_hat;
            elseif alpha_correlation <= eps_correlation && beta_correlation > eps_correlation
                rho = beta_hat;
            else
                %rho = rho; %do nothing
            end
            end
            
            %update gamma
            if adjust_gamma == 1
            if alpha_correlation>eps_correlation && beta_correlation>eps_correlation
                gamma=1+(2*sqrt(alpha_hat*beta_hat))/(alpha_hat+beta_hat);
            elseif alpha_correlation>eps_correlation && beta_correlation<=eps_correlation
                gamma=1.9;
            elseif alpha_correlation<=eps_correlation && beta_correlation>eps_correlation
                gamma=1.1;
            else
                gamma=1.5;
            end
            end
            
            l_hat_0 = l_hat;
            l_0     = l;
            z_0     = z;
            x_0     = x;
        end
        
end

end

%% Sopping conditions
function [stop,adjust_rho,ind_ref]=stop_admm(log,i,evol_rel_tol,feas_tol,obj_tol,adjust_rho,ind_ref)
    stop=0;
    if i>25 && log.set_feasibility(end)<feas_tol && abs(log.obj(end)-log.obj(end-1))/log.obj(end-1) < obj_tol
        disp(['stationary objective and reached feasibility, exiting ADMM (iteration ',num2str(i),' )'])
        stop=1;
    end
    if i>1 && log.evol_x(i)<evol_rel_tol
        disp(['relative evolution to small, exiting ADMM (iteration ',num2str(i),' )'])
        stop=1;
    end
    if i>1 && log.r_pri(i)>max(log.r_pri((i-1):-1:max((i-50),1))) && adjust_rho==0
        disp(['no primal residual reduction, fixing rho (iteration ',num2str(i),' )'])
         adjust_rho = 0;
         ind_ref = i;
    elseif adjust_rho == 0 && i>(ind_ref+10) && log.r_pri(i)>max(log.r_pri((i-1):-1:max(ind_ref,max((i-50),1)))) 
       disp(['no primal residual reduction, exiting ADMM (iteration ',num2str(i),' )'])
       stop = 1;
    end 
end
    
%% Output checks
function []=output_check_ADMM(Ax,funProj,x,i)

%feasibility check
nAx      = norm(Ax);
feas_err = norm(funProj(Ax)-Ax)/nAx;

    disp(['iteration: ',num2str(i),',  final feasibility error ADMM : ',num2str(feas_err)])


if isreal(x)==0
    disp('warning: Result of ADMM is not real')
end
if nnz(x)~=length(x)
    disp('warning: Result of ADMM contains ZEROS')
end

end