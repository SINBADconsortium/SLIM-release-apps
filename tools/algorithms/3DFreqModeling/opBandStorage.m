classdef opBandStorage < opSpot & handle
%OPBANDSTORAGE Band-storage matrix. Multiply and divide operators supported.
%  Since this operator uses iterative methods such as CARPCG,CARPCR,etc. for solving linear
%     systems, and a subset of these methods need to scale the rows of the matrix, the user
%     has the option of what options this SPOT operator should support (forward / adjoint modes,
%     multiplication / division operations) so that the minimal amount of information to perform
%     them is stored. Typically, if you are interested in solving the helmholtz equation in an FWI
%     context, you are rarely looking to apply the multiplication mode, so the corresponding
%     matrix entries for forward/adjoint multiplication are irrelevant.
%
%
% Curt Da Silva, 2015
%
% Usage:
%
%   A = opBandStorage(coef,idx,opts);
%
% Input:
%   coef              - bandwidth x N matrix of coefficients (short and fat matrix)
%   idx               - 1 x bandwidth matrix of coefficient offsets
%   opts              - parameter struct
%     .mode           - store the minimum amount of information to operate in one of the following modes
%                       opBandStorage.FORW_MULT_ONLY     - forward multiplication only
%                       opBandStorage.FORW_DIV_ONLY      - forward division only
%                       opBandStorage.MULT_ONLY          - forward/adjoint multiplication only
%                       opBandStorage.DIV_ONLY           - forward/adjoint division only
%                       opBandStorage.FORW_MULT_DIV      - forward multiplication/division only
%                       opBandStorage.MULT_DIV           - both multiplication and division operations supported
%
% if division is chosen
%     .solve_opts     - a LinSolveOpts object, specifying the linear
%                       solver parameters for the operator
%
% Output:
%   A                 - SPOT operator
    
    properties (SetAccess = private)
        mode,
        coef, idx,
        coef_adj,idx_adj,    
        row_norms_sq,
        adj_row_norms_sq,
        mvp_nthreads;
    end
    properties 
        solve_opts,top_level;
    end
    properties (Constant)
        FORW_MULT_ONLY = 'forw_mult';
        FORW_DIV_ONLY  = 'forw_div';        
        MULT_ONLY      = 'forw_adj_mult';
        DIV_ONLY       = 'forw_adj_div';
        FORW_MULT_DIV  = 'forw_mult_div';
        MULT_DIV       = 'mult_div';
    end
    methods(Static)
        function p = is_mult_mode(mode)
            % True if the specified string corresponds to a
            % multiplication-supported mode of opBandStorage, false
            % otherwise
            p = opBandStorage.is_valid_op(mode,'forw_mult') || opBandStorage.is_valid_op(mode,'adj_mult');
        end
        
        function p = is_div_mode(mode)
            % True if the specified string corresponds to a
            % division-supported mode of opBandStorage, false
            % otherwise
            p = opBandStorage.is_valid_op(mode,'forw_div') || opBandStorage.is_valid_op(mode,'adj_div');
        end
        
        function p = is_valid_op(mode,operation)
            switch operation
              case 'forw_mult'
                p = length(find(ismember(mode,{opBandStorage.FORW_MULT_ONLY,opBandStorage.MULT_ONLY, opBandStorage.FORW_MULT_DIV, opBandStorage.MULT_DIV }))) > 0;
              case 'adj_mult'
                p = length(find(ismember(mode,{opBandStorage.MULT_ONLY, opBandStorage.MULT_DIV }))) > 0;
              case 'forw_div'
                p = length(find(ismember(mode,{opBandStorage.FORW_DIV_ONLY, opBandStorage.DIV_ONLY, opBandStorage.FORW_MULT_DIV, opBandStorage.MULT_DIV }))) > 0;
              case 'adj_div'
                p = length(find(ismember(mode,{opBandStorage.DIV_ONLY, opBandStorage.MULT_DIV }))) > 0;
              case 'adjoint'
                p = length(find(ismember(mode,{opBandStorage.MULT_ONLY,opBandStorage.DIV_ONLY, opBandStorage.MULT_DIV }))) > 0;

              otherwise 
                error('Unrecognized operation');
            end
        end        
    end
    methods
        function op = opBandStorage(coef,idx,opts)
            Nt = size(coef,2);
            op = op@opSpot('Band-storage matrix',Nt,Nt);
            
            assert(exist('opts','var')~=0, 'Need to specify options');
            
            op.coef = []; op.idx = []; op.coef_adj = []; op.idx_adj = [];
            opts.mode = opBandStorage.MULT_ONLY;
            op.sweepflag = true;
            op.top_level = false;
            op.cflag = ~isreal(coef);
            op.mvp_nthreads = 1;
            is_in = @(x,opts) ~isempty(find(ismember(x,opts)));
            
            % Determine what we need to keep and what we can ignore
            forw_mult = is_in(opts.mode,{opBandStorage.FORW_MULT_ONLY,opBandStorage.MULT_ONLY,opBandStorage.FORW_MULT_DIV,opBandStorage.MULT_DIV});
            adj_mult  = is_in(opts.mode,{opBandStorage.MULT_ONLY,opBandStorage.MULT_DIV});
            forw_div  = is_in(opts.mode,{opBandStorage.FORW_DIV_ONLY,opBandStorage.DIV_ONLY,opBandStorage.FORW_MULT_DIV,opBandStorage.MULT_DIV});
            adj_div   = is_in(opts.mode,{opBandStorage.DIV_ONLY,opBandStorage.MULT_DIV});
            op.mode = opts.mode;

            if forw_mult
                op.coef = coef; op.idx = idx;
            end
            
            if adj_mult
                [op.coef_adj,op.idx_adj] = Htransp(coef,idx);
                op.coef_adj = conj(op.coef_adj);
            end             
            
            if forw_div || adj_div
                assert(isa(opts.solve_opts,'LinSolveOpts'),'solve_opts must be a LinSolveOpts object');
                op.solve_opts = copy(opts.solve_opts);
                if strcmp(op.solve_opts.solver,LinSolveOpts.SOLVE_CGMN)||strcmp(op.solve_opts.solver,LinSolveOpts.SOLVE_CRMN)
                    op.solve_opts.precond = LinSolveOpts.PREC_KACZSWP;
                end
                if ~forw_mult && forw_div
                    op.coef = coef; op.idx = idx;
                end
                if ~adj_mult && adj_div
                     [op.coef_adj,op.idx_adj] = Htransp(coef,idx);
                     op.coef_adj = conj(op.coef_adj);
                end
                
                switch op.solve_opts.precond
                  case {LinSolveOpts.PREC_KACZSWP,LinSolveOpts.PREC_MLCR}                   
                    op.row_norms_sq = vec(sum(abs(coef).^2,1));
                    if adj_div
                        op.adj_row_norms_sq = vec(sum(abs(op.coef_adj).^2,1));
                    else
                        op.adj_row_norms_sq = [];
                    end
                  otherwise
                    op.row_norms_sq = []; op.adj_row_norms_sq = [];
                end
            else
                op.solve_opts = [];
            end                                                  
        end
        
        function x = kacz_sweep(op,w,x,b,nsweeps,mode)
            % KACZ_SWEEP - Kaczmarz sweeping preconditioner for solving Ax = b
            %
            % Usage:
            %   y = A.kacz_sweep(w,x,b,nsweeps,R,idx);
            %
            % Input:
            %   A         - opBandStorage matrix
            %   w         - over relaxation parameter
            %   x         - current model vector
            %   b         - right hand side of Ax = b
            %   nsweeps   - number of sweeps to perform
            %   R, idx    - band storage coefficients + indices to use
            %
            % Output:
            %   y         - kaczmarz sweep applied to x
            %            
                if mode==1
                    R = op.coef; idx = op.idx; row_norms_sq = op.row_norms_sq;
                else
                    R = op.coef_adj; idx = op.idx_adj; row_norms_sq = op.adj_row_norms_sq;
                end
                for k=1:nsweeps
                    if true
                        x = sweep_MT_mex(R,idx,x,b,w,row_norms_sq);
                    else
                        % matlab implementation of the kaczmarz sweep, for reference but not actual use
                        N = size(R,2); 
                        for m=1:2*N
                            if m > N, i = 2*N-m+1; else i = m; end
                            I = i+idx; j = find(I >=1 & I <= N);
                            I = I(j);
                            c = w * ( b(i,:) - (R(j,i)).'*x(I,:))./row_norms_sq(i);
                            x(I,:) = x(I,:) + conj(R(j,i))*c;
                        end
                    end
                end            
            end
        
        function set_solve_opts(op,solveopts)
            assert(isa(solveopts,'LinSolveOpts'));
            op.solve_opts = solveopts;
        end
        
        function y = linsolve(op,x,mode,x0,maxit,solver,precond)
            if mode == 1 && ~opBandStorage.is_valid_op(op.mode,'forw_div')
                error('Forward division unsupported');
            elseif mode ~= 1 && ~opBandStorage.is_valid_op(op.mode,'adj_div')
                error('Adjoint division unsupported');
            end
            I = find(sqrt(sum(abs(x).^2,1)) > 1e-12);
            y = zeros(size(x));
            if isempty(I), return; end
            if length(I) < size(x,2)
                x = x(:,I); x0 = x0(:,I);            
            end
            lsopts = copy(op.solve_opts);
            if exist('maxit','var')
                lsopts.maxit = maxit;
            end
            if exist('solver','var')
                lsopts.solver = solver;
            end
            %Parse + setup preconditioner
            if strcmp(lsopts.solver,LinSolveOpts.SOLVE_CRMN)||strcmp(lsopts.solver,LinSolveOpts.SOLVE_CGMN)
                lsopts.precond = LinSolveOpts.PREC_KACZSWP;       
            elseif exist('precond','var')
                lsopts.precond = precond;                                                
            end
            
            if ~isa(lsopts.precond,'function_handle') && ~isa(lsopts.precond,'opSpot')
                switch lsopts.precond
                  case LinSolveOpts.PREC_KACZSWP
                    w = 1.5; nsweeps = 1;
                    zb = 0*x;                                         
                    b = kacz_sweep(op,w,zb,x,nsweeps,mode); 
                    A = opDirac(op.n)-opFunction_swp(op.n,op.n,@(y,m) kacz_sweep(op,w,y,zb,nsweeps,mode),op.cflag);
                    lsopts.precond = @(x) x;
                  otherwise
                    lsopts.precond = @(x) x;
                    A = op; b = x;
                end
            else          
                b = x; 
                if mode == 1
                    A = op; 
                else
                    A = op'; 
                    if isa(lsopts.precond,'opSpot')
                        lsopts.precond = lsopts.precond';
                    end
                end
            end
            
            % Setup solver
            switch lsopts.solver
              case LinSolveOpts.SOLVE_CGMN
                solve = @(A,b) CGMN(A,b,x0,lsopts);              
              case LinSolveOpts.SOLVE_CRMN                                                   
                solve = @(A,b) CRMN(A,b,x0,lsopts);
              case LinSolveOpts.SOLVE_FGMRES     
                %for debugging purposes
                if op.top_level
                    if ~isprop(lsopts,'top_level'), lsopts.addprop('top_level'); end
                    lsopts.top_level = true;
                    lsopts.addprop('output_freq');
                    lsopts.output_freq = 1;
                end

                if size(x,2) == 1
                    solve = @(A,b) FGMRES(A,b,x0,lsopts);
                else
                    solve = @(A,b) FGMRES_MT(A,b,x0,lsopts);
                end
              otherwise
                error('Unknown solver flag');
            end
            
            out = solve(A,b);
            if length(I) < size(b,2);
                y(:,I) = out;
            else
                y = out;
            end
           
        end
    end
    
    methods (Access = protected)
        function y = multiply(op,x,mode)
            if mode==1
                if ~opBandStorage.is_valid_op(op.mode,'forw_mult'), error('Invalid mode set - forward multiplication not supported'); end
                if true                    
                    y = Hmvp(op.coef,op.idx,x,op.mvp_nthreads);
                else
                    % matlab implementation of band-storage
                    % matrix-vector product
                    y = zeros(size(x));
                    N = size(op,2);
                    for m=1:N
                        i = m;                    
                        I = i+op.idx; j = find(I >=1 & I <= N);
                        I = I(j);
                        y(i) = op.coef(j,i).'*x(I,:);
                    end  
                end
            else
                if ~opBandStorage.is_valid_op(op.mode,'adj_mult'), error('Invalid mode set - adjoint multiplication not supported'); end
                y = Hmvp(op.coef_adj,op.idx_adj,x,op.mvp_nthreads);
            end
            
        end
        
        function y = divide(op,x,mode)
            y = linsolve(op,x,mode,zeros(size(x)));
        end
    end
    
end

