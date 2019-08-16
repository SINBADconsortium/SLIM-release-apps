% Multi Level CRMN due to R. Lago and F. Herrmann (to appear)
% 
% This is a new preconditioner to be used by FGMRES. It is inspired by T2V 
% due to Calandra et. al, but uses one iteration of CRMN as smoother on
% every level, and 5 iterations of CRMN as coarse solver.
%
% This routine merely sets it up and returns  a pointer to the proper function,
% as it requires quite a lot of parameters setting and it might be a bit too
% much to expect the user to set them all manually; to be used by (flexible)
% GMRES,  BiCGStab or any non-Hermitian Krylov solver.
% 
% The preconditioner is composed of a three grid V cycle, but it requires 4
% operators - two Helmholtz and two shifted Laplacian:
% 
%                           I_h     I_h
%    (V cycle on             \      /     
%     Helmholtz operator)     \    /
%                              \  /          
%                              I_2h  ==>  I_2h   I_2h       
%                                          \     /     
%                                           \   /      (V cycle on shifted 
%                                            \ /        Laplacian operator) 
%                                           I_4h       
%
% In this version, all the operators are computed based on the finest one, that
% needs to be provided as the structure helm_fine.
% 
% IMPORTANT:
% Unlike GMG_CGPV12_T2V, this preconditioner is supposed to work equally well 
% for 7 or 27 points stencils.
% 
% Usage:  
%  precond = ML_CRMN(Hk,v,comp_grid,model,freq,opts,flog)
%            
% 
% INPUT:
%   Hk         - discretized Helmholtz operator (SPOT operator or matrix)
%   v          - current model corresponding to Hk
%   comp_grid  - computational grid structure (output of discretize_helmholtz)
%   model      - the same structure you passed to discrete_helmholtz routine
%                in order to obtain a discrete operator.
%   freq       - 
%   opts       - 
%
%
%   flog       - [OPTIONAL] file identifier for the (already open) log file. 
%                Set to 0 (or simply do not provide it) if no verbosity is 
%                desired.
%
% OUTPUT:
%   precond     - a pointer to the preconditioner function. 
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: March, 2015
%
% Updated: Curt Da Silva, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function precond = ML_CRMN(Hk,v,comp_grid,model,freq,opts,flog)

if ~exist('flog','var')
   flog = 0;
end

% Preconditioner parameters
% level 1 - finest level
lv1_smoother                           = LinSolveOpts.SOLVE_CRMN;
lv1_smoother_prec                      = LinSolveOpts.PREC_KACZSWP;
lv1_smoother_maxit                     = 1;
lv1_coarse_solver_params.maxit         = 1; % # outer iterations of FGMRES
lv1_coarse_solver_params.maxinnerit    = 5; % # inner iterations of FGMRES

% number of V cycles for finest level
par_lv1.maxit = 1;

% level 2 - semi coarse level
lv2_smoother                           = LinSolveOpts.SOLVE_CRMN;
lv2_smoother_prec                      = LinSolveOpts.PREC_KACZSWP;
lv2_smoother_maxit                     = 1;
lv2_coarse_solver                      = LinSolveOpts.SOLVE_CRMN;
lv2_coarse_solver_prec                 = LinSolveOpts.PREC_IDENTITY;

lv2_coarse_maxit                       = 20;

% number of level 2 V cycles
par_lv2.maxit                          = 1;

interp_basis = 'cubic';
%Shifted Laplacian constants
beta1 =  1;
%beta2 = .5*sqrt(-1);
beta2 = 0;
adj_mode = opBandStorage.is_valid_op(opts.mode,'adjoint');

% Opts struct for semi coarsest level
opts_lv2         = opts;
opts_lv2.dt      = 2*comp_grid.dt;
opts_lv2.pml     = ceil(comp_grid.pml.x/2);
opts_lv2.scheme  = comp_grid.scheme;
if isa(opts.solve_opts,'matlab.mixin.Copyable')
    opts_lv2.solve_opts  = copy(opts.solve_opts);    
end
opts_lv2.solve_opts.precond = LinSolveOpts.PREC_KACZSWP;
opts_lv2.solve_opts.solver = LinSolveOpts.SOLVE_CRMN;

% Opts struct for coarsest level
opts_lv3  = opts;
opts_lv3.dt      = 4*comp_grid.dt;
opts_lv3.pml     = ceil(comp_grid.pml.x/4);
opts_lv3.scheme  = comp_grid.scheme;
if isa(opts.solve_opts,'matlab.mixin.Copyable')
    opts_lv3.solve_opts  = copy(opts.solve_opts);    
end
opts_lv3.solve_opts.precond = LinSolveOpts.PREC_KACZSWP;
opts_lv3.solve_opts.solver = LinSolveOpts.SOLVE_CRMN;

% Coarsen model
[lv1_to_lv2,~,n_lv2] = fine2coarse(comp_grid.nt_nopml,comp_grid.dt,opts_lv2.dt,interp_basis);
v_lv2 = lv1_to_lv2*vec(v);

[lv2_to_lv3] = fine2coarse(n_lv2,opts_lv2.dt,opts_lv3.dt,interp_basis);
v_lv3 = lv2_to_lv3*vec(v_lv2);


%% Build up Level 2 (semi coarsest) MG preconditioner
% Level 2 shifted Laplacian
[SLap_lv2,comp_grid_lv2] = discrete_helmholtz(v_lv2,model,sqrt(beta1+beta2)*freq,opts_lv2,flog); 

% Level 3 shifted Laplacian (coarsest)
[SLap_lv3,comp_grid_lv3] = discrete_helmholtz(v_lv3,model,sqrt(beta1+beta2)*freq,opts_lv3,flog); 

[lv2_to_lv3,lv3_to_lv2] = fine2coarse(comp_grid_lv2.nt,comp_grid_lv3.nt,interp_basis);

par_lv2.restriction  = @(u)(lv2_to_lv3*u); 
par_lv2.prolongation = @(u)(lv3_to_lv2*u);

par_lv2.pre_smoother = @(b,x0)(linearsolve(SLap_lv2,b,1,x0,lv2_smoother_maxit,lv2_smoother,lv2_smoother_prec));
par_lv2.pst_smoother = par_lv2.pre_smoother;
par_lv2.coarse       = @(b,x0)(linearsolve(SLap_lv3,b,1,x0,lv2_coarse_maxit,lv2_coarse_solver,lv2_coarse_solver_prec));
level2_precon        = @(b)(GMG_V(SLap_lv2,b,0*b,par_lv2));

lvl2_precon = opFunction_swp(size(SLap_lv2,2),size(SLap_lv2,2),{level2_precon,@(x) x},1);

if adj_mode
    % copy solver settings for the adjoint level 2 solver
   par_lv2_adj = struct; par_lv2_adj.restriction = par_lv2.restriction; par_lv2_adj.prolongation = par_lv2.prolongation;
   par_lv2_adj.maxit = par_lv2.maxit;
   par_lv2_adj.pre_smoother = @(b,x0) linearsolve(SLap_lv2,b,-1,x0,lv2_smoother_maxit,lv2_smoother,lv2_smoother_prec); par_lv2_adj.pst_smoother = par_lv2_adj.pre_smoother;
   par_lv2_adj.coarse = @(b,x0) linearsolve(SLap_lv3,b,-1,x0,lv2_coarse_maxit,lv2_coarse_solver,lv2_coarse_solver_prec);
   level2_adj_precon = @(b) GMG_V(SLap_lv2',b,0*b,par_lv2_adj);
   
end


%% Build up level 1 preconditioner (finest)

[lv1_to_lv2,lv2_to_lv1] = fine2coarse(comp_grid.nt,comp_grid_lv2.nt,interp_basis);

par_lv1.restriction  = @(u)(lv1_to_lv2*u); 
par_lv1.prolongation = @(u)(lv2_to_lv1*u);

lv1_coarse_solver_params.tol       = -1;
lv1_coarse_solver_params.precond   = lvl2_precon;

% Assemble everything - Level 1
%-----------------------------

par_lv1.pre_smoother = @(b,x0)(linearsolve(Hk,b,1,x0,lv1_smoother_maxit,lv1_smoother,lv1_smoother_prec));

Helm_lv2 = discrete_helmholtz(v_lv2,model,freq,opts_lv2,flog);
par_lv1.coarse       = @(b,x0)(FGMRES(Helm_lv2,b,x0,lv1_coarse_solver_params));
par_lv1.pst_smoother = par_lv1.pre_smoother;

prec = @(b)(GMG_V(Hk,b,[],par_lv1));

if adj_mode
    % copy solver settings for the adjoint mode - lvl 1 solver
    lv1_coarse_solver_params_adj.tol = -1; lv1_coarse_solver_params_adj.precond = level2_adj_precon;
    lv1_coarse_solver_params_adj.maxit = lv1_coarse_solver_params.maxit; 
    lv1_coarse_solver_params_adj.maxinnerit = lv1_coarse_solver_params.maxinnerit;    
    par_lv1_adj.restriction = par_lv1.restriction; par_lv1_adj.prolongation = par_lv1.prolongation;
    par_lv1_adj.pre_smoother = @(b,x0) linearsolve(Hk,b,-1,x0,lv1_smoother_maxit,lv1_smoother,lv1_smoother_prec);
    
    par_lv1_adj.coarse = @(b,x0) FGMRES(Helm_lv2',b,x0,lv1_coarse_solver_params_adj);
    par_lv1_adj.pst_smoother = par_lv1_adj.pre_smoother;
    par_lv1_adj.maxit = par_lv1.maxit;
    
    prec_adj = @(b) GMG_V(Hk',b,[],par_lv1_adj);
    
    precond = opFunction_swp(size(Hk,1),size(Hk,1),{prec,prec_adj},1);
else
    precond = opFunction_swp(size(Hk,1),size(Hk,1),{prec,@(x) x},1);
end


end
