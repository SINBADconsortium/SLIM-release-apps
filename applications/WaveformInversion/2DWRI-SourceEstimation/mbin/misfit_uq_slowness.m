function [f g H w] = misfit_uq_slowness(m,D,Q,model,params);

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

[fmis gmis Hmis w] = fh(m(:));


[fpr gpr Hpr]    = fhp(m(:));

f    = fmis + beta * fpr;
g    = gmis + beta * gpr;
H    = Hmis + beta * Hpr;

if model.Priormp > 0
    if isfield(model,'velprior')
        velprior = model.velprior;
    else
        velprior = 0;
    end
    if velprior == 0
        mp       = model.mp;
        SIGMA_mp = model.SIGMA_mp;
        fhpm     = @(x) misfit_prior_mp(x,mp,SIGMA_mp);
        [fpm gpm Hpm] = fhpm(m(:));
    else
        mp_v       = model.mp;
        %mp_v     = 1./sqrt(mp);
        SIGMA_mp = model.SIGMA_mp;
        fhpm     = @(x) misfit_prior_mp(x,mp_v,SIGMA_mp);
        dfv      = @(x) -0.5 * x(:).^(-1.5);
        Jv       = opDiag(dfv(m(:)));
        
        [fpm gpm Hpm] = fhpm(1./sqrt(m(:)));
        gpm           = Jv' * gpm;
        Hpm           = Jv' * Hpm * Jv + opDiag(gpm) * opDiag(3/4 * m(:).^(-2.5));
    end
    f        = f + fpm;
    g        = g + gpm;
    H        = H + Hpm;
    
end









