classdef opMultigrid < opSpot & handle
% opMultigrid - Multigrid operator
%
    properties (SetAccess = protected)
        nlevels, smoother, coarse_solver;
        H,smoothers,restriction,prolongation,solve_coarse;
        shape;
    end
    
    properties (Constant)
        ML_GMRES = 'ml_gmres';
    end
    
    properties
        recursive_v_cycle;
    end
    
    methods
        function op = opMultigrid(nlevels,shape)
            if exist('shape','var')==0 || isempty(shape)
                shape = 'V';
            end
            op = op@opSpot('Multigrid operator',0,0);
            op.nlevels = nlevels;            
            op.smoother = [];
            op.coarse_solver = [];            
            op.H = [];
            op.smoothers = [];
            op.restriction = [];
            op.prolongation = [];
            op.solve_coarse = [];
            op.recursive_v_cycle = false;
            op.shape = shape;
        end
        function op = set_smoother(op,smoother)
            if isa(smoother,'LinSolveOpts')
                op.smoother = smoother;
            else
                error('smoother must be an object of type LinSolveOpts');
            end
        end
        function op = set_coarse_solver(op,solver)
            if isa(solver,'LinSolveOpts')
                op.coarse_solver = solver;
            else
                error('solver must be an object of type LinSolveOpts');
            end
        end
        function op1 = extract_subv_cycle(op,level)
            op1 = opMultigrid(op.nlevels-(level-1));
            op1.set_smoother(op.smoother);
            op1.set_coarse_solver(op.coarse_solver);
            op1.H = op.H(level:end);
            op1.m = size(op1.H{1},1); op1.n = size(op1.H{1},2);
            op1.smoothers = op.smoothers(level:end);
            op1.restriction = op.restriction(level:end);
            op1.prolongation = op.prolongation(level:end);
            op1.solve_coarse = op.solve_coarse;
            
        end
        function op = construct_external(op,Hs,S,R,P,coarse_solve)
            op.n = size(Hs{1},1);
            op.m = op.n;
            op.H = Hs;
            op.smoothers = S;
            op.restriction = R;
            op.prolongation = P;
            op.solve_coarse = coarse_solve;
        end
        
        function op = construct(op,H,v,comp_grid,model,freq,opts,recursive_v_cycle)
        % discretize helmholtz on different levels
            if isempty(op.smoother) || isempty(op.coarse_solver)
                error('Smoother, coarse solver must be set before calling this method');
            end
            opts.solve_opts = copy(opts.solve_opts);
            opts.disp_output = false;
            
            % dummy values, won't be used
            opts.solve_opts.solver = LinSolveOpts.SOLVE_FGMRES;
            opts.solve_opts.precond = LinSolveOpts.PREC_IDENTITY;            
            
            % grid size parameters
            
            if exist('recursive_v_cycle','var')==0
                recursive_v_cycle = false;
            end

            op.recursive_v_cycle = recursive_v_cycle;
            op.n = size(H,2); op.m = size(H,1);
            
            [H,smoother,restriction,prolongation,coarse_solver] = construct_multigrid(H,v,comp_grid,model,freq,opts,op.smoother,op.coarse_solver,op.nlevels);
            op.H = H; op.smoothers = smoother; op.restriction = restriction; op.prolongation = prolongation; op.solve_coarse = coarse_solver;
            if recursive_v_cycle
                gmg_opts = struct; gmg_opts.maxit = 1;
                ls_opts = copy(op.coarse_solver);                 
                coarse_solve = LinearSolver(op.H{end},ls_opts);
                for i=op.nlevels-1:-1:1
                    Hf = op.H{i}; Hc = op.H{i+1};
                    gmg_opts.restriction = op.restriction{i}; 
                    gmg_opts.prolongation = op.prolongation{i};
                    gmg_opts.pre_smoother = op.smoothers{i};
                    ls_opts = copy(op.coarse_solver);
                    ls_opts.precond = coarse_solve;
                    gmg_opts.coarse = LinearSolver(Hc,ls_opts);
                    coarse_solve = FuncObj(@GMG_V,{Hf,[],[],gmg_opts,[]},{'A','b','x','opts','mode'});
                end
                op.solve_coarse = coarse_solve;
            end
        end
    end
    
    methods (Access = protected)
        function x = multiply(op,b,mode)
            x = zeros(size(b));
            if op.recursive_v_cycle      
                x = op.solve_coarse(b,0*b,mode);
            else            
                x_lvl = cell(1,op.nlevels+1);
                x_lvl{1} = x;
                b_lvl = cell(1,op.nlevels+1);
                b_lvl{1} = b;
                % recursive smoothing, restriction
                for i=1:op.nlevels-1                    
                    S = op.smoothers{i};
                    x = S(b_lvl{i},x,mode);
                    if mode==1
                        r = b_lvl{i} - op.H{i}*x;
                    else
                        r = b_lvl{i} - op.H{i}'*x;
                    end
                    R = op.restriction{i}; 
                    x_lvl{i+1} = R*r;
                    x = x_lvl{i+1};
                    b_lvl{i+1} = R*b_lvl{i};
                end
                % coarse solve
                x = op.solve_coarse(x,0*x,mode);
                
                switch op.shape 
                  case 'V'                                        
                    % Do nothing
                    
                  case 'F'
                    for j=1:op.nlevels-1
                        
                    end
                    x = op.solve_coarse(x,0*x,mode);
                  otherwise
                    error(['Unrecognized shape ' op.shape]);
                end
                
                % recursive prolongation, smoothing
                for i=op.nlevels-1:-1:1                
                    Pr = op.prolongation{i};
                    x = x_lvl{i} + Pr*x;
                    S = op.smoothers{i};
                    x = S(b_lvl{i},x,mode);
                end    
                    
            end
        end
    end
    
end