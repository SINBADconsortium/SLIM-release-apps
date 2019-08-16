function [w,taxis] = wvlet(f,dt,t);
% ricker wavelet difination
% input:
% 	f: central frequency
% 	dt: time interval
% 	t: total sampling time
%
% Author: Xiang Li
	
nw=6./f/dt;
nw=2*floor(nw/2)+1;
nc=floor(nw/2);
i=1:nw;
alpha=(nc-i+1).*f*dt*pi;
beta=alpha.^2;
w=(1.-beta.*2).*exp(-beta);


nt = round(t./dt)+1;
if length(w) < nt
	w(end:nt) = 0;
end

taxis = (0:nt-1) * dt;
