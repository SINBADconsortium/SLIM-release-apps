function [f g H H2] = misfit_prior_smth_vel(x, sigma_p, n, delta);
    %% Function to calculate the misfit for the non smooth prior information of the 
    %  Bayesian inversion
    %  Usage: 
    %       [f g H] = misfit_prior_smth_vel(x, sigma_p, n);
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
    
    [m dm ddm]       = slowsquare2vel(x(:));
    Hdm              = opDiag(dm);
    
    if nargin < 4
        [fvel gvel Hvel Hvel2] = misfit_prior_smth(m, 1, n);
   else
        [fvel gvel Hvel Hvel2] = misfit_prior_smth(m, 1, n, delta);
   end
    
    
    
    lambda = 1/sigma_p^2;
    dvel   = gvel .* ddm(:);
    
    f    = lambda * fvel;
    g    = lambda * (gvel .* dm(:));
    H    = lambda * (Hdm * Hvel * Hdm + opDiag(dvel));
    
    H2   = sqrt(lambda) * [Hvel2*Hdm; opDiag(sqrt(dvel))];
    
    
    
    
    
    
