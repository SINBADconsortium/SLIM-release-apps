function w = getrick(f,o,dt,nt)
% w = getrick(f,o,dt,nt)
% w = (1-2*(pi*f*(t-o))^2)*exp(-(pi*(t-o))^2);
t = (0:nt-1)*dt;
r = (pi*f*(t-o));
w = (1-2*r.^2).*exp(-r.^2);
