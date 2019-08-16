function op = opFFTsym(n)
% opFFTsym  One-dimensional symmetric real fast Fourier transform (FFT).
%
%    opFFTsym(N) create a one-dimensional normalized Fourier transform along the 
%    operator for vectors of length N.
%
%    The output will only contain the positive frequencies (as well as DC and possibly Nyquist).
%

%   Copyright 2008, Tim Lin

% calculate output vector length
if mod(n,2) == 1
    m = ceil(n/2);
    has_nyquist = 0;
else
    m = (n/2) + 1;
    has_nyquist = 1;
end

subfunc_handle = @(x,mode) opFFTsym_intrnl(x,mode);

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction


function y = opFFTsym_intrnl(x,mode)

    % actual operation
    if mode == 0
       y = {m,n,[0,1,0,1],{'FFTsym'}};
    elseif mode == 1
       y = fft(x) / sqrt(n);
       y(m+1:end) = [];
       % doubling the magnitude of everything except DC and Nyquist, because the conjugate freqs are chopped off
       y(2:end-1) = 2*y(2:end-1);
       if not(has_nyquist), y(end) = 2*y(end); end
    else
       y = ifft(x,n,'symmetric') * sqrt(n);
    end
end

end