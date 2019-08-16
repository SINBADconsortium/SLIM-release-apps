function val = SNR(xref,x)
%SNR - SIGNAL TO NOISE RATIO
%
% Usage:
%   val = SNR(xref,x);
%
% Input:
%   xref - reference signal
%   x    - test signal
% 
% Output:
%   val  - SNR
%

val = -20*log10(norm(vec(xref) - vec(x))/norm(vec(xref)));

end