function w = getrickRomb(f,o,dt,nt,Rtemp)

% getrickPhaseEncod(model.f0,model.t0,dt,nt,model.phase(:,ls))
% w = getrick(f,o,dt,nt)
% w = (1-2*(pi*f*(t-o))^2)*exp(-(pi*(t-o))^2);
t = (0:nt-1)*dt;
r = (pi*f*(t-o));
tw = (1-2*r.^2).*exp(-r.^2);

w  = Rtemp*tw';
