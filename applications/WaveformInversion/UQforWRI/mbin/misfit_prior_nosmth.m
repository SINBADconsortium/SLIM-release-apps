function [f g H] = misfit_prior_nosmth(x, sigma_p, n);
    %% Function to calculate the misfit for the non smooth prior information of the 
    %  Bayesian inversion
    %  Usage: 
    %       [f g H] = misfit_prior_nosmth(x, sigma_p, n);
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
    
    