function [Hmis Lpr Lpm] = misfit_uq_slowness_smp(m,D,Q,model,params);
    
    %% function of the misfit for the uncertainty quantification
    % Use : [f g H] = misfit_uq_slowness(m,D,Q,model,params);
    % Input : 
    % m - model 
    % D - Data
    % Q - Source
    % model - model parameters
    % params - computing parameters
    % 
    % Output:
    % f  - misfit function
    % g  - gradient
    % H  - Hessian
    % 
    % Author : Zhilong Fang
    % Date   : 2016/01
    
    sigma           = model.sigma;
    beta            = model.beta;
    sigma_p         = model.sigma_p;
    for i = 1:length(model.freq)
        params.lambda(i) = params.lambda(i) * sigma(i);
    end    
    
    fh      = misfit_setup(Q,D,model,params);
    
    if ~isfield(model,'flagPriorVel')
        flagPriorVel = 0;
    else
        flagPriorVel = model.flagPriorVel;
    end
    
    switch model.PriorType
        case 'nonsmooth' 
            if flagPriorVel < 1
                fhp     = @(x) misfit_prior_nosmth(x, sigma_p, model.n);
            else
                fhp     = @(x) misfit_prior_nosmth_vel(x, sigma_p, model.n);
            end
        case 'smooth'
            if flagPriorVel < 1
                fhp     = @(x) misfit_prior_smth(x, sigma_p, model.n);
            else
                fhp     = @(x) misfit_prior_smth_vel(x, sigma_p, model.n);
            end
        case 'TV'
            if flagPriorVel < 1
                fhp     = @(x) misfit_prior_TV(x, sigma_p, model.n);
            else
                fhp     = @(x) misfit_prior_TV_vel(x, sigma_p, model.n);
            end
        otherwise
            error('Wrong model.PriorType, please select one of the following types: nonsmooth, smooth, TV');
    end
    
    [fmis gmis Hmis] = fh(m(:));
    
    
    [fpr gpr Hpr Lpr]    = fhp(m(:));
    
    f    = fmis + beta * fpr;
    g    = gmis + beta * gpr;
    H    = Hmis + beta * Hpr;
    
    if model.Priormp > 0
        mp       = model.mp;
        SIGMA_mp = model.SIGMA_mp;
        fhpm     = @(x) misfit_prior_mp(x,mp,SIGMA_mp);
        [fpm gpm Hpm Lpm] = fhpm(m(:));
        f        = f + fpm;
        g        = g + gpm;
        H        = H + Hpm;
    end
    

    
    
    
    
    
    
    
    