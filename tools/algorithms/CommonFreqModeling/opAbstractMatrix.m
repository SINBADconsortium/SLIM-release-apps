classdef opAbstractMatrix < opSpot & handle
% Abstract matrix, supports explicit or implicit matrix-vector multiplication. For explicit
% matrices that need to be inverted, computes + stores LU factors. For implicit matrices, 
% stores FuncObj objects to produce forward/adjoint multiplications, jacobi sweeps, etc.
% 
% Curt Da Silva, 2016
%
% Usage:
%   A = opAbstractMatrix(matrix_mode,params,solve_opts);
%
% Input:
%   matrix_mode - one of 'explicit' or 'implicit'
%   params      - struct of 
%    in explicit mode:
%         .coef   - matrix coefficients
%    
%    in implicit mode:
%         .N         - size of the matrix, N x N
%         .multiply  - FuncObj or function_handle with signature params.multiply(x,mode);
%         .jacobi    - FuncObj or function_handle with signature params.jacobi(alpha,b,x,mode);
%         .kacz_sweep- FuncObj or function_handle with signature params.kaczswp(w,x,b,nsweeps,mode);
%         .iscomplex - whether the underlying matrix is complex or not
%  
%   solve_opts   - LinSolveOpts object, for solving linear systems, or empty if only multiply mode desired
    
    properties ( SetAccess = protected )
        mode,params;
    end
    properties
        solve_opts;
    end
    methods
        function op = opAbstractMatrix(matrix_mode,params,solve_opts)                        
            switch matrix_mode
              case 'explicit'                
                Nt = size(params.coef,2);
                iscomplex = ~isreal(params.coef);
                is_full = ~issparse(params.coef);
                params.mode = 'explicit';
                params.is_full = is_full;
                if strcmp(solve_opts.solver,LinSolveOpts.SOLVE_LU)
                    if is_full
                        [L,U,P] = lu(params.coef);
                        params.L = L; params.U = U; params.P = P;
                    else
                        [L,U,P,Q,R] = lu(params.coef);
                        params.L = L; 
                        params.U = U;
                        params.P = P; 
                        params.Q = Q;
                        params.R = R;
                    end
                end
              case 'implicit'
                assert(isfield(params,'multiply') && (isa(params.multiply,'FuncObj')||isa(params.multiply,'function_handle')),'params.multiply must be of type FuncObj or a function_handle');
                assert(isfield(params,'N'), 'params.N must be present');
                assert(isfield(params,'iscomplex'),'params.iscomplex must be present');
                Nt = params.N;
                iscomplex = params.iscomplex;
              otherwise
                error('Unrecognized matrix format');
            end
            op = op@opSpot('Abstract Matrix',Nt,Nt);
            op.cflag = iscomplex;
            op.sweepflag = false;
            
            op.solve_opts = solve_opts;                        
            op.mode = matrix_mode;
            op.params = params;
        end        
        function y = jacobi(op,alpha,b,x,mode,maxiter)
        % Jacobi smoother iteration
            switch op.mode
              case 'explicit'                
                y = jacobi(op.params.coef,alpha,b,x,mode,maxiter);
              case 'implicit'
                if isfield(op.params,'jacobi') && ~isempty(op.params.jacobi)
                    for i=1:maxiter
                        x = op.params.jacobi(alpha,b,x,mode);
                    end
                    y = x;
                else
                   error('No jacobi method specified for this matrix');
                end
            end
        end
        function y = kaczswp(op,w,x,b,nsweeps,mode)
        % Kaczmarz sweep
            switch op.mode
              case 'explicit'
                y = kaczswp(op.params.coef,w,x,b,nsweeps,mode);
              case 'implicit'
                if isfield(op.params,'kacz_sweep') && ~isempty(op.params.kaczswp)
                    y = op.params.kacz_sweep(w,x,b,nsweeps,mode);
                else
                    error('No kaczmarz sweep specified for this matrix');
                end                
            end
        end
    end
    methods (Access = protected)
        function y = multiply(op,x,mode)
            switch op.mode
              case 'explicit'
                if mode == 1
                    y = op.params.coef*x;
                else
                    y = op.params.coef'*x;
                end
              case 'implicit'
                y = op.params.multiply(x,mode);
            end
        end
        function y = divide(op,x,mode)
            if norm(x) < 1e-12, y= zeros(size(x));
            else
                if isempty(op.solve_opts) || ~isa(op.solve_opts,'LinSolveOpts')
                    error('op.solve_opts is empty or an invalid class');
                end
                y = linearsolve(op,x,[],mode,op.solve_opts);
            end
        end
    end
end