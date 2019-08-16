function [f g H H2] = misfit_prior_smth(x, sigma_p, n);
    %% Function to calculate the misfit for the non smooth prior information of the 
    %  Bayesian inversion
    %  Usage: 
    %       [f g H] = misfit_prior_smth(x, sigma_p, n);
    %  Input:
    %  x          - vectorized 2D model
    %  sigma_p    - standard deviation for the nonsmooth distribution
    %  n          - size of the model
    %  
    %  Output:
    %  f          - misfit function value
    %  g          - gradient
    %  H          - Hessian
    %
    %  Author:
    %  Zhilong Fang, SLIM, UBC
    %  2016/01
    
    lambda = 1/sigma_p^2;
    
    az0  = ones(n(1),1);
    az1  = [-.5*ones(n(1)-2,1);-1;-1];
    az2  = [-1;-1;-.5*ones(n(1)-2,1)];
    az   = [az1 az0 az2]/2;
    Dz   = spdiags(az, [-1 0 1], n(1),n(1));
    Ix   = speye(n(2));
    
    ax0  = ones(n(2),1);
    ax1  = [-.5*ones(n(2)-2,1);-1;-1];
    ax2  = [-1;-1;-.5*ones(n(2)-2,1)];
    ax   = [ax1 ax0 ax2]/2;
    Dx   = spdiags(ax, [-1 0 1], n(2),n(2));
    Iz   = speye(n(1));
    
    S1   = opKron(Ix, Dz);
    S2   = opKron(Dx, Iz);
    S    = (S1 + S2);
    
    y    = S * x;
    
    f    = .5*lambda * norm(y)^2;
    g    = lambda  * S' * y;
    H    = lambda  * S' * S;
    H2   = sqrt(lambda) * S;
    
    
    
    
    
    
