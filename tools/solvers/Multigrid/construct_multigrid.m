function [Hs,smoothers,restriction,prolongation,coarse_solver] = construct_multigrid(H,v,comp_grid,model,freq,opts,ls_smoother,ls_coarse,nlevels)
% CONSTRUCT_MULTIGRID - creates all of the necessary multigrid components for inverting the helmholtz matrix
%
% Curt Da Silva, 2015
% 
% Usage:
%   [Hs,smoothers,restriction,prolongation,coarse_solver] = construct_multigrid(H,v,comp_grid,model,freq,opts,ls_smoother,ls_coarse,nlevels);
%
%
    coarse_factor = 2;
    opts.solve_opts = copy(opts.solve_opts);
    % dummy values, won't be used
    opts.solve_opts.solver = LinSolveOpts.SOLVE_FGMRES;
    opts.solve_opts.precond = LinSolveOpts.PREC_IDENTITY;            
            
    % grid size parameters
    nt_nopml_fine = comp_grid.nt_nopml;
    nt_fine = comp_grid.nt;
    dt_fine = comp_grid.dt;
    
    pml_fine = comp_grid.pml(1,1);
    mode_2d = length(nt_fine)==2 || nt_fine(3)==1;

    Hs = cell(nlevels,1);
    Hs{1} = H;
    smoothers = cell(nlevels-1,1);
    restriction = cell(nlevels-1,1);
    prolongation = cell(nlevels-1,1);
    
    for i=1:nlevels-1
        opts.dt = coarse_factor^i*dt_fine;

        opts.pml = ceil(pml_fine/coarse_factor^i);
        if mode_2d, opts.dt(3) = 1; end
        % Smoothers                                
        if isa(ls_smoother,'cell') 
            smoothers{i} = LinearSolver(H,ls_smoother{i});
        else
            smoothers{i} = LinearSolver(H,ls_smoother);
        end
        % Coarsen model
        [to_coarse,~,n_coarse] = fine2coarse(nt_nopml_fine,opts.dt/coarse_factor,opts.dt);
        if mode_2d, assert(length(n_coarse)==2 || n_coarse(3)==1); end
        v_coarse = reshape(to_coarse*vec(v),n_coarse);
        if i==nlevels-1,opts.solve_opts = ls_coarse; end
        [H,comp_grid_coarse] = discrete_helmholtz(v_coarse,model,freq,opts);
        
        % Interpolation operators with PML nodes
        [to_coarse,to_fine,nt_coarse] = fine2coarse(nt_fine,comp_grid_coarse.nt,'linear');        
        
        % weighted restriction
        restriction{i} = 1/8*to_coarse;
        
        % linear interpolation
        prolongation{i} = to_fine;
        
        Hs{i+1} = H;
                
        v = v_coarse; nt_nopml_fine = n_coarse; nt_fine = nt_coarse;
    end
    
    coarse_solver = LinearSolver(Hs{end},ls_coarse);
end
