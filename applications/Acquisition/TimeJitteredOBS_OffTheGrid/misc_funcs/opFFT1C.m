function op = opFFT1C(n)

% OPFFT1C  One-dimensional fast Fourier transform (FFT).
%
% Copyright 2008, Gilles Hennenfent
%
% July, 2012 : Wrapped into a SPOT operator
%              Haneet Wason


fh = @(x, mode) opFFT1C_intrnl(n, x, mode);
op = opFunction(n, n, fh);


function y = opFFT1C_intrnl(n, x, mode)
if mode == 0
    y = {n, n, [1,1,1,1], {'FFT1C'}};
elseif mode == 1
    % Synthesis 
    y = fftshift(ifft(ifftshift(x))) * sqrt(length(x));
else
    % Analysis
    y = fftshift(fft(ifftshift(x))) / sqrt(length(x));
end

