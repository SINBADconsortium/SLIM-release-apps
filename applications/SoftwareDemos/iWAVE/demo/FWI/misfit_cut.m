function [f g] = misfit_cut(x,model,data,options,ncut);

if nargout < 2
     f = misfit(x,model,data(:),options);
else
     [f g] = misfit(x,model,data(:),options);
end
f           = f*1e6;
g1       = reshape(g(1:prod(model.n)),model.n);
g1(1:ncut,:) = 0;
g(1:prod(model.n)) = g1;
g = g*1e6;
