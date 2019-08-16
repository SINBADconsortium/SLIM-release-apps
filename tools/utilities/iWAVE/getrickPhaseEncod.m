function w = getrickPhaseEncod(f,o,dt,nt,teta_loc)

% getrickPhaseEncod(model.f0,model.t0,dt,nt,model.phase(:,ls))
% w = (1-2*(pi*f*(t-o))^2)*exp(-(pi*(t-o))^2);
t  = (0:nt-1)*dt;
r  = (pi*f*(t-o));
tw = (1-2*r.^2).*exp(-r.^2);

%w = sign(real(exp(i*teta_loc(:))))'.*real(ifft(fft(tw).*exp(i*teta_loc')));
%w = sign(real(exp(i*teta_loc(1))))*real(ifft(fft(tw).*exp(i*teta_loc')));

exptemp = fft(real(ifft(exp(i*teta_loc'))));
exptemp = exptemp./abs(exptemp);
w       = sign(real(exp(i*teta_loc(1))))*real(ifft(fft(tw).*exptemp));
