function Mtop = ML_CALANDRA12(H,v,comp_grid,model,freq,opts)
% Multi-level Helmholtz preconditioned by the Shifted Laplacian, implementation of Algorithm 3 in 
% 'An improved two-grid preconditioner for the solution of three-dimensional Helmholtz problems in heterogeneous media'
% H. Calandra, S. Gratton, X. Pinel, X. Vasseur, NUMERICAL LINEAR ALGEBRA WITH APPLICATIONS, 2013
%
% Curt Da Silva, 2016
% 
% Usage:
%   M = ML_CALANDRA12(H,v,comp_grid,model,freq,opts);
%
%
% 
    nlevels = 4;
    jacobi_s = cell(nlevels,1);
    
    % Set up smoothers for each level
    % Parameters given by equation (16)
    alphas = [0.8, 0.8, 0.2, 1];
    for i=1:nlevels
        jacobi_s{i} = LinSolveOpts();
        jacobi_s{i}.solver = LinSolveOpts.SOLVE_FGMRES;
        jacobi_s{i}.maxit = 1;
        jacobi_s{i}.maxinnerit = 10;
    end
    
    % Coarse solver for Shifted Laplacian
    coarse_s = LinSolveOpts();
    coarse_s.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_s.maxit = 1;
    coarse_s.maxinnerit = 10;
    coarse_s.tol = 5e-1;
    coarse_s.precond = jacobi_s{nlevels};
    
    beta = 0.5;
    
    [SL,S,R,P,C] = construct_multigrid(H,v,comp_grid,model,sqrt(1+beta*1i)*freq,opts,jacobi_s,coarse_s,nlevels);
    
    %Shifted laplacian preconditioner for the 2h level
    M = opMultigrid(nlevels-1);
    construct_external(M,SL(2:end),S(2:end),R(2:end),P(2:end),C);
    
    % Parameters in Section 5.1
    smoother = LinSolveOpts();
    smoother.solver = LinSolveOpts.SOLVE_FGMRES;
    smoother.maxit = 3;
    smoother.maxinnerit = 5;
    %smoother.precond = jacobi_s{1};
    
    coarse_s = LinSolveOpts();
    coarse_s.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_s.maxit = 1;
    coarse_s.maxinnerit = 20;
    coarse_s.precond = M;
    coarse_s.tol = 5e-1;
    
    %Top level multigrid operator
    Mtop = opMultigrid(2);
    set_smoother(Mtop,smoother);
    set_coarse_solver(Mtop,coarse_s);
    construct(Mtop,H,v,comp_grid,model,freq,opts,false);
    
end

function Mtop = ML_CALANDRA123(H,v,comp_grid,model,freq,opts)
% Multi-level Helmholtz preconditioned by the Shifted Laplacian, implementation of Algorithm 3 in 
% 'An improved two-grid preconditioner for the solution of three-dimensional Helmholtz problems in heterogeneous media'
% H. Calandra, S. Gratton, X. Pinel, X. Vasseur, NUMERICAL LINEAR ALGEBRA WITH APPLICATIONS, 2013
%
% Curt Da Silva, 2016
% 
% Usage:
%   M = ML_CALANDRA12(H,v,comp_grid,model,freq,opts);
%
%
% 
    nlevels = 4;
    jacobi_s = cell(nlevels,1);
    
    % Set up smoothers for each level
    % Parameters given by equation (16)
    alphas = [0.8, 0.8, 0.2, 1];
    for i=1:nlevels
        jacobi_s{i} = LinSolveOpts();
        jacobi_s{i}.solver = LinSolveOpts.SOLVE_JACOBI;
        jacobi_s{i}.maxit = 2;
        jacobi_s{i}.params.alpha = alphas(i);
    end
    
    % Coarse solver for Shifted Laplacian
    coarse_s = LinSolveOpts();
    coarse_s.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_s.maxit = 1;
    coarse_s.maxinnerit = 10;
    coarse_s.precond = jacobi_s{nlevels};
    
    beta = -0.5;
    
    [SL,S,R,P,C] = construct_multigrid(H,v,comp_grid,model,sqrt(1+beta*1i)*freq,opts,jacobi_s,coarse_s,nlevels);
    
    %Shifted laplacian preconditioner for the 2h level
    M = opMultigrid(nlevels-1);
    construct_external(M,SL(2:end),S(2:end),R(2:end),P(2:end),C);
    
    % Parameters in Section 5.1
    smoother = LinSolveOpts();
    smoother.solver = LinSolveOpts.SOLVE_FGMRES;
    smoother.maxit = 1;
    smoother.maxinnerit = 2;
    %smoother.precond = jacobi_s{1};
    
    coarse_s = LinSolveOpts();
    coarse_s.solver = LinSolveOpts.SOLVE_FGMRES;
    coarse_s.maxit = 2;
    coarse_s.maxinnerit = 10;
    coarse_s.precond = M;
    coarse_s.tol = 5e-1;
    
    %Top level multigrid operator
    Mtop = opMultigrid(2);
    set_smoother(Mtop,smoother);
    set_coarse_solver(Mtop,coarse_s);
    construct(Mtop,H,v,comp_grid,model,freq,opts,false);
    
end
