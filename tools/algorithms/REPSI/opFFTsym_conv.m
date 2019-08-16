function op = opFFTsym_conv(n)
% opFFTsym_conv  One-dimensional symmetric real fast Fourier transform (FFT) that does not include self-adjoint scaling
%                and does not double the magnitude of the positive frequencies. WILL NOT PASS DOT-TEST.
%                ONLY meant for computing convolutions where the FFT appears with a corresponding inverse.
%                In this case the adjoint (mode = 2) actually implements the inverse.
%
%    opFFTsym_conv(N) create a one-dimensional normalized Fourier transform along the 
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

subfunc_handle = @(x,mode) opFFTsym_conv_intrnl(x,mode);

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction


function y = opFFTsym_conv_intrnl(x,mode)

    % actual operation
    if mode == 0
       y = {m,n,[0,0,0,0],{'FFTsym_conv'}};
    elseif mode == 1
       y = fft(x);
       y(m+1:end) = [];
    else
       y = ifft(x,n,'symmetric');
    end
end

end

