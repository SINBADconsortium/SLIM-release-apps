% Three level Multigrid Preconditioner for 2D wavefield solves
%
%  
%
% Usage:  
%  T = ML_GMRES2D(H,v,comp_grid,model,freq,opts)
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
function T = ML_GMRES2D(H,v,comp_grid,model,freq,opts)
    
    nlevels = 3;
    T = opMultigrid(nlevels);
    smoother = LinSolveOpts();
    smoother.solver = LinSolveOpts.SOLVE_FGMRES;
    smoother.precond = LinSolveOpts.PREC_IDENTITY;
    smoother.maxit = 1; smoother.maxinnerit = 10;
    smoother.tol = 1e-6;
    
    coarse_solver = LinSolveOpts();
    coarse_solver.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_solver.precond = LinSolveOpts.PREC_IDENTITY;
    coarse_solver.maxit = 1; 
    coarse_solver.maxinnerit = 10;
    coarse_solver.tol = 5e-1;
    
    set_smoother(T,smoother);
    set_coarse_solver(T,coarse_solver);
    recursive_v_cycle = true;
    construct(T,H,v,comp_grid,model,freq,opts,recursive_v_cycle);
   
end
