function [b,f] = fftrl(a,t,mode)
% fft for real-valued vectors
%
% Tristan van Leeuwen, 2012
% tleeuwen@eos.ubc.ca
%
% use:
%   [b,f] = fftrl(a,t,mode)
%
% input:
%   a - input data
%   t - time vector
%   mode: 1: forward, -1:inverse
%
% output:
%   b - vector

nt = length(t);
dt = t(2) - t(1);
if mod(nt,2)==0
	nf=nt/2;
else
	nf = floor(nt/2) + 1;
end
tmax = t(end) - t(1);
f = 0:1/tmax:.5/dt;

switch mode
    case 1
        b = fft(a,[],1);
        b = b(1:nf,:);
    case -1
        a = [a;conj(a(ceil(nt/2):-1:2,:))];
        b = ifft(a);
        b = real(b);
    otherwise
        error('Unknown mode');
end
