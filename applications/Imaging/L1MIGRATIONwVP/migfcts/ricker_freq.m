function w = ricker_freq(f0,t0,freq)
% compute spectrum of a ricker wavelet
% f0: peak freq; t0: time shift; freq: freq. samples
% by Tuning, SLIM, 2014
    
w = exp(1i*2*pi*freq*t0);
w = (freq).^2.*exp(-(freq/f0).^2).*w;
    