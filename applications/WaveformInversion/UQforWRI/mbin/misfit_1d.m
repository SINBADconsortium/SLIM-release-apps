function [f g H] = misfit_1d(m,Q,D,model,params,Px)

sigma           = model.sigma;
for i = 1:length(model.freq)
    params.lambda(i) = params.lambda(i) * sigma(i);
end

fh      = misfit_setup(Q,D,model,params);

[fmis gmis Hmis] = fh(Px*m(:));

mp       = model.mp;
SIGMA_mp = model.SIGMA_mp;
fhpm     = @(x) misfit_prior_mp(x,mp,SIGMA_mp);
[fpm gpm Hpm] = fhpm(m(:));

f        = fmis + fpm;
g        = Px'*gmis + gpm;
H        = Px'*Hmis*Px + Hpm;
