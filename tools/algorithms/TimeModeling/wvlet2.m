function [w,taxis] = wvlet2(f,dt,t,t0)
% function [w,taxis] = wvlet(f,dt,t,t0)
% ricker wavelet difination
% input:
% 	f: central frequency
% 	dt: time interval
% 	t: total sampling time
%
% Author: Xiang Li


nw     =6./f/dt;
nw     =2*floor(nw/2)+1;
nc     =floor(nw/2);
i      =1:nw;
alpha  =(nc-i+1).*f*dt*pi;
beta   =alpha.^2;
w_temp =(1.-beta.*2).*exp(-beta);

nt    = round(t./dt)+1;

if nargin < 4
    idx = find(abs(w_temp)>eps);
    t0  = (length(w_temp)./2 - idx(1)).*dt;
end

t_temp = ((0:(length(w_temp)-1))-floor(length(w_temp)./2)).* dt + t0;

taxis = (0:nt-1) * dt;

Pt     = opLInterp1D(t_temp,taxis);

w      = Pt * w_temp(:);
