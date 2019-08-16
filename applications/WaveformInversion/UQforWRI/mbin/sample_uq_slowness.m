function [Smp] = sample_uq_slowness(m,D,Q,model,params);
    
    %% function of using Gibbs sampling method to sample
    % Use : [Smp] = misfit_uq_slowness(m,D,Q,model,params);
    % Input : 
    % m - model 
    % D - Data
    % Q - Source
    % model - model parameters
    % params - computing parameters
    % 
    % Output:
    % Smp  - samplers
    % 
    % Author : Zhilong Fang
    % Date   : 2016/01
    
    sigma           = model.sigma;
    beta            = model.beta;
    sigma_p         = model.sigma_p;
    nsmp            = params.nsmp;
    Smp             = zeros(prod(model.n), nsmp);
    np              = prod(model.n);
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
    
    if model.Priormp > 0
        mp       = model.mp;
        SIGMA_mp = model.SIGMA_mp;
        fhpm     = @(x) misfit_prior_mp(x,mp,SIGMA_mp);
    end
    
    for s_idx = 1:nsmp
    
        [fmis gmis Hmis] = fh(m(:));
    
    
        [fpr gpr Hpr]    = fhp(m(:));
        
        Hpr              = double(Hpr);
        Hpr              = sparse(Hpr);
    
        f    = fmis + beta * fpr;
        g    = gmis + beta * gpr;
        H    = Hmis + beta * Hpr;
    
        if model.Priormp > 0
            [fpm gpm Hpm] = fhpm(m(:));
            Hpm           = double(Hpm);
            Hpm           = sparse(Hpm);
        end
        f        = f + fpm;
        g        = g + gpm;
        H        = H + Hpm;
        HL       = chol(H);
        m        = m - HL\(HL'\g) + HL'\randn(np,1);
        
        Smp(:,s_idx) = m(:);
        s_idx
        % if mod(s_idx,10) == 1;
        %     figure;imagesc(reshape(m,model.n))
        %     keyboard
        % end
        
    end

    
    
    
    
    
    
    
    