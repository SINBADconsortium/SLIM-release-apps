classdef LinSolveOpts2D < LinSolveOpts
%LINSOLVEOPTS2D Default parameters for solving 2D helmholtz
%equations (using Matlab's sparse LU decomposition)
%
% Author: Curt Da Silva
%  
    
    methods
        function opts = LinSolveOpts2D()
            opts.solver     = LinSolveOpts.SOLVE_LU;
        end
    end
    
end

