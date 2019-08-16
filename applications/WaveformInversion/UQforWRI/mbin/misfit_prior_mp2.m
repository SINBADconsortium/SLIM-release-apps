function [f g H L] = misfit_prior_mp(m,mp,ISIGMA)
    %% function to calculate the prior misfit with a prior model mp
    % Usage:
    %       [f g H] = misfit_prior_mp(m,mp,SIGMA)
    % Input:
    % m  - current model
    % mp - prior model
    % SIGMA - standard deviation vector or matrix
    % 
    % Output
    % f  - misfit value
    % g  - gradient
    % H  - Hessian
    %
    % Author : Zhilong Fang, SLIM, UBC
    % Date   : 2016/01
    
   
    dm     = m(:) - mp(:);
    f      = 0.5 * dm' * (ISIGMA* dm);
    g      = ISIGMA* dm;
    H      = ISIGMA;
%    L      = opDiag(sqrt(lambda));
