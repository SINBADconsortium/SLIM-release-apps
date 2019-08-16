function [stdx, cix, Smp] = RTO_Sampling(H, H2, LAll, betaAll, gammaAll, mmap, params);
    %% function to use Randomized then optimize sampling method
    % Usage:
    % [stdx, civ, Smp] = RTO_Sampling(H, H2, LAll, sigmaAll, mmap, params);
    %
    % Input:
    % H       - Hessian of the misfit function
    % H2      - square root of the Hessian
    % LAll    - Structure that contains all prior information matrix
    % betaALL - Structure that contains all penalty parameters for 
    %           the prior information
    % gammaAll- Structure that contains all relative model for 
    %           the prior information
    %           beta * (m - gamma)^{T} * L'L * (m - gamma)
    % mmap    - maximum a posterior point,
    % params  - necessary parameters
    %         - nsmp        : number of samplers
    %         - alpha_CI    : Confidence level
    %         - tol         : tolerance for pcg
    %         - maxiter     : maximal iterations for pcg
    %
    % Outpur
    % stdx    - standard deviation of samplers
    % cix     - confidence interval of samplers
    % Smp     - all generated samplers
    %
    %  Author:
    %  Zhilong Fang, SLIM, UBC
    %  2016/01 
    
    n        = size(H,1);
    nsmp     = params.nsmp;
    alpha_CI = params.alpha_CI;
    tol      = params.tol;
    maxiter  = params.maxiter;
    [nH1 nH2] = size(H2);
    
    
    ymmap    = H * mmap;
    Hfull    = H;
    
    for i = 1: length(betaAll)
        Li      = LAll{i};
        betai   = betaAll{i};
        gammai  = gammaAll{i};
        Hfull   = Hfull + betai * Li' * Li;
    end
    
    for smp_j = 1:nsmp
        y   = ymmap + H2' * randn(nH1,1);
        for i = 1: length(betaAll)
            Li      = LAll{i};
            betai   = betaAll{i};
            gammai  = gammaAll{i};
            si      = sqrt(betai) * Li * gammai + randn(size(Li,1),1);
            yt      = sqrt(betai) * real(Li' * si);
            y       = y + yt;           
        end
        Smp(:,smp_j) = pcg(Hfull, y, tol, maxiter);
        if isfield(params,'tmp_name');
            tmp_name = params.tmp_name;
        else
            tmp_name = 'TMP_SAMPLE';
        end
        save(tmp_name,'Smp','-v7.3'); 
        if mod(smp_j,5) == 1
            fprintf(1,'%6.2d\n',smp_j);
        end
        
    end
    
    for i = 1: size(Smp,1);
        stdx(i)  = std(Smp(i,:));
        cix(i,:) = confidence_interval(Smp(i,:),alpha_CI);
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
