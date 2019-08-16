function [f g H H2] = misfit_prior_smth(x, sigma_p, n, delta);
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
    
    if nargin < 4
        delta = [x(1)*0.01; x(end)*0.01];
    end
    
    lambda = 1/sigma_p^2;
    
    az0  = ones(n(1),1);
    az1  = [-.5*ones(n(1)-2,1);-1;-1];
    az2  = [-1;-1;-.5*ones(n(1)-2,1)];
    az   = [az1 az0 az2]/2;
    Dz   = spdiags(az, [-1 0 1], n(1),n(1));
    
    Dz(1,2) = 0;
    Dz(1,1) = 1/delta(1) * sigma_p;
    Dz(end,end-1) = 0;
    Dz(end,end)     = 1/delta(2) * sigma_p;     
    
    Ix   = speye(n(2));
    
    ax0  = ones(n(2),1);
    ax1  = [-.5*ones(n(2)-2,1);-1;-1];
    ax2  = [-1;-1;-.5*ones(n(2)-2,1)];
    ax   = [ax1 ax0 ax2]/2;
    Dx   = spdiags(ax, [-1 0 1], n(2),n(2));
    Dx(1,2) = 0;
    Dx(1,1) = 1/delta(1) * sigma_p;
    Dx(end,end-1) = 0;
    Dx(end,end)     = 1/delta(2) * sigma_p;     
    Iz   = speye(n(1));
    
    S1   = opKron(Ix, Dz);
    SS1 = kron(Ix, Dz);
    S2   = opKron(Dx, Iz);
    SS2 = kron(Dx, Iz);
    S    = (S1 + S2);
    SS = SS1 + SS2;
    
    y    = S * x;
    
    f    = .5*lambda * norm(y)^2;
    g    = lambda  * S' * y;
    H    = lambda  * S' * S;
    H2   = sqrt(lambda) * SS;
    
    
    
    
    
    
