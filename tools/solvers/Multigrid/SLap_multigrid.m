function M = SLap_multigrid(v,comp_grid,model,freq,opts)
% Multigrid operator for the Shifted Laplacian preconditioner
% 
% Curt Da Silva, 2016
%
% Usage:
%   M = SLap_multigrid(v,comp_grid,model,freq,opts);
%
% Input:
%   v         - velocity model
%   comp_grid - computational grid struct
%   model     - model struct
%   freq      - frequency
%   opts      - opts struct
%  
% Output:
%   M         - opMultigrid, M*q approximately inverts the Shifted Laplacian corresponding to H
%
    nlevels = 4;
    jacobi_s = cell(nlevels,1);
    alphas = [0.8, 0.8, 0.2, 1];
    for i=1:nlevels
        jacobi_s{i} = LinSolveOpts();
        jacobi_s{i}.solver = LinSolveOpts.SOLVE_JACOBI;
        jacobi_s{i}.maxit = 2;
        jacobi_s{i}.params.alpha = alphas(i);
    end
    
    coarse_s = LinSolveOpts();
    coarse_s.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_s.maxit = 1;
    coarse_s.maxinnerit = 10;
    coarse_s.precond = jacobi_s{nlevels};    

    beta = -0.5;
    opts_loc = opts;
    opts_loc.solve_opts = copy(opts.solve_opts);
    opts_loc.solve_opts.solver = LinSolveOpts.SOLVE_FGMRES;
    opts_loc.solve_opts.precond = LinSolveOpts.PREC_IDENTITY;
    
    SLtop = discrete_helmholtz(v,model,(sqrt(1+beta*1i))*freq,opts_loc);
                
    [SL,S,R,P,C] = construct_multigrid(SLtop,v,comp_grid,model,sqrt(1+beta*1i)*freq,opts,jacobi_s,coarse_s,nlevels);
    
    M = opMultigrid(nlevels);
    construct_external(M,SL,S,R,P,C);
    
end
