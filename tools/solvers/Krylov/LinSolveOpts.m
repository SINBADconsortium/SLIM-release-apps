classdef LinSolveOpts < dynamicprops & matlab.mixin.Copyable
%LINSOLVEOPTS Parameter object for solving linear systems
%
% Author: Curt Da Silva
%  
% Options
%   .tol          - relative residual tolerance (default: 1e-6)
%
%   .maxit        - maximum (outer) iterations (default: 10000)
%
%   .maxinnerit   - maximum (inner) iterations, for certain solvers (default: 10)
%
%   .solver       - linear solver to use, one of
%      LinSolveOpts.CGMN    - CGMN 
%      LinSolveOpts.CRMN    - CRMN
%      LinSolveOpts.GMRES   - GMRES
%      LinSolveOpts.FGMRES  - GMRES w/ a flexible preconditioning option
% 
%   .precond      - preconditioner to use, one of
%      LinSolveOpts.PREC_KACZSWP  - kaczmarz sweeps (default for CGMN, CRMN)
%      LinSolveOpts.PREC_IDENTITY - identity preconditioner (default)
%      LinSolveOpts.PREC_MLCR     - multi-level CR preconditioner    
%      LinSolveOpts.PREC_MLGMRES  - multi-level GMRES preconditioner
%  
%      OR 
%      a LinSolveOpts object, which specifies an iterative solver to use 
%     
    
    properties
        tol, maxit, maxinnerit, solver, precond,params;
    end
    
    properties (Constant)
        SOLVE_CGMN     = 'solve_cgmn';
        SOLVE_CRMN     = 'solve_crmn';
        SOLVE_GMRES    = 'solve_gmres';        
        SOLVE_FGMRES   = 'solve_fgmres';
        SOLVE_BICGSTAB = 'solve_bicgstab';
        SOLVE_LU       = 'solve_lu';
        SOLVE_BACKSLASH= 'solve_backslash';
        SOLVE_JACOBI   = 'solve_jacobi';
        PREC_KACZSWP   = 'prec_kaczsweep';
        PREC_MLCR      = 'prec_mlcr';
        PREC_MLCG      = 'prec_mlcg';
        PREC_T2V       = 'prec_t2v';
        PREC_IDENTITY  = 'prec_identity';        
        PREC_BICGSTAB  = 'prec_bicgstab';
        PREC_MLGMRES   = 'prec_mlgmres';
        PREC_SHIFTLAP  = 'prec_shiftlap';
        PREC_CALANDRA12 = 'prec_calandra12';
    end
    
    methods
        function opts = LinSolveOpts()
            opts.tol        = 1e-6;
            opts.maxit      = 10000;
            opts.maxinnerit = 5;
            opts.solver     = LinSolveOpts.SOLVE_FGMRES;
            opts.precond    = LinSolveOpts.PREC_IDENTITY;
            opts.params     = struct;
        end
    end
    
end

