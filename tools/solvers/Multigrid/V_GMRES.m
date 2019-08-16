% Three level Multigrid Preconditioner
%
%  
%
% Usage:  
%  T = ML_GMRES(H,v,comp_grid,model,freq,opts)
%            
% 
% INPUT:
%   H          - discretized Helmholtz operator (SPOT operator or matrix)
%   v          - current model corresponding to H
%   comp_grid  - computational grid structure (output of discretize_helmholtz)
%   model      - the same structure you passed to discrete_helmholtz routine
%                in order to obtain a discrete operator.
%   freq       - frequency (Hz)
%   opts       - opts struct, passed in by discrete_helmholtz
%
%
% OUTPUT:
%   T          - ML GMRES preconditioner in a SPOT operator (T corresponds to H, T' corresponds to H')
%
% AUTHOR: Curt Da Silva, 2015
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function T = V_GMRES(H,v,comp_grid,model,freq,opts)
    ks_outer = 1;
    ks_inner = 5;
    kc_outer = 1;
    kc_inner = 5;
    coarsetol = 5e-1;
    
    if isfield(opts,'k')
        ks_outer = opts.k(1);
        ks_inner = opts.k(2);
        kc_outer = opts.k(3);
        kc_inner = opts.k(4);
    end
    
    
    
    nlevels = 2;
    T = opMultigrid(nlevels);
    smoother = LinSolveOpts();
    smoother.solver = LinSolveOpts.SOLVE_FGMRES;
    smoother.precond = LinSolveOpts.PREC_IDENTITY;
    smoother.maxit = ks_outer; 
    smoother.maxinnerit = ks_inner;
    smoother.tol = 1e-6;
    
    coarse_solver = LinSolveOpts();
    coarse_solver.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_solver.precond = LinSolveOpts.PREC_IDENTITY;
    coarse_solver.maxit = kc_outer; 
    coarse_solver.maxinnerit = kc_inner;
    coarse_solver.tol = coarsetol;
    
    set_smoother(T,smoother);
    set_coarse_solver(T,coarse_solver);
    recursive_v_cycle = false;
    construct(T,H,v,comp_grid,model,freq,opts,recursive_v_cycle);
    
end
