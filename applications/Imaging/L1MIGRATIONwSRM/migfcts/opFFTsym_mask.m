function op = opFFTsym_mask(n, mask)
% Syntax:
% op = opFFTsym_mask(n, mask)
% "mask" is optional
%
% Description:
% One-dimensional symmetric real fast Fourier transform (FFT), 
% supporting frequency mask. The output will only contain the 
% positive frequencies (as well as DC and possibly Nyquist).
%
% Input list:
% n: length of the vector
% mask: (optional input) frequency mask, MUST be logical if present
%
% Output list:
% op: the operator
%
% Author: Ning Tu (Adapted from Tim's opFFTsym)
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% calculate output vector length
if mod(n,2) == 1
    m = ceil(n/2);
    has_nyquist = 0;
else
    m = (n/2) + 1;
    has_nyquist = 1;
end

if not(exist('mask','var'))
    mask = true(m,1);
end
if not(islogical(mask))
    mask = logical(mask);
    disp(['Warning: mask is not of logical type. Forced logical ' ...
          'conversion.'])
end
compressed_size = sum(double(mask));
if compressed_size > m
    error('Fatal: frequency mask length larger than Nyquist');
end

subfunc_handle = @(x,mode) opFFTsym_mask_intrnl(x,mode);
op = opFunction(compressed_size,n,subfunc_handle);

function y = opFFTsym_mask_intrnl(x,mode)

% actual operation
    if mode == 0
        y = {compressed_size,n,[0,1,0,1],{'FFTsym_mask'}};
    elseif mode == 1
        y = fft(x) / sqrt(n);
        y(m+1:end) = [];
        % doubling the magnitude of everything except DC and
        % Nyquist, because the conjugate freqs are chopped off
        y(2:end-1) = 2*y(2:end-1);
        if not(has_nyquist)
            y(end) = 2*y(end);
        end
        y = y(mask);
    else
        y = zeros(m,1);
        y(mask) = x;
        y = ifft(y,n,'symmetric') * sqrt(n);
    end
end

end