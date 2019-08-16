function [f g] = funMoreauYosida(r, params)
rho = params.rho;
pro = max(0, abs(r) - rho).*sign(r);
f = norm(pro, 1) + 0.5*norm(r - pro)^2/rho;
g = (r-pro)./rho;
end