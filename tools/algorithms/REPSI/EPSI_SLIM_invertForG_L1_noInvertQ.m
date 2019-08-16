function [x, residual, residual_list, iters_performed, sigma_reached] = EPSI_SLIM_invertForG_L1_onInvertQ(A, b, sigma, xprev, tau, options)
% L1 inversion for primary impulse response (using SPGL1) when not simultaneously inverting for Q (only change iteration definitions)
    
    SPGL1_EXIT_ROOT_FOUND = 1; % const
    
    spgl1_opts  = spgSetParms('verbosity' ,     options.verbosity*2, ...
                        'bpTol'     ,     1e-5, ...
                        'decTol'    ,     options.decTol, ...
                        'minPareto',      options.minL1iters, ...
                        'ignorePErr',     options.spgl1_ignorePErr, ...
                        'subspaceMin',    options.subspaceMin, ...
                        'iteration',      options.maxTotalIter, ...
                        'nPrevVals' ,     4);
                        
    [x, r, g, info] = spgl1(A, b, tau, sigma, xprev, spgl1_opts);
    clear r
    clear g
    
    % final residual ebergy
    residual = info.rNorm;
    % per-iteration residual energy
    residual_list = info.rNorm2(:);
    % number of iterations executed
    iters_performed = info.iter;
    % check if mismatch between data and prediction wavefields are below acceptable levels
    sigma_reached = (info.stat == SPGL1_EXIT_ROOT_FOUND);
    
    % exact line search scaling of the solution
    g = A' * b;
    y = A * x;    
    s = (x' * g) / (y' * y);
    
    s = abs(undist(s)); disp(['   ... sacled by: ' num2str(s)])
    x = s .* x;
    
end