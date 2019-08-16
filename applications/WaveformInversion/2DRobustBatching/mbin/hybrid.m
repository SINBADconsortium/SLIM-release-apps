function [f g h] = hybrid(x, m)
% Hybrid misfit: sqrt(1 + x.^2/m)
%
% use:
%   [f,g,h] = hybrid(x,{m})
%
% input:
%   x - vector  
%   m - scale
%
% output:
%   f - huber misfit
%   g - gradient
%   h - hessian
%

% Author: Tristan van Leeuwen, Aleksandr Aravkin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% input checking
if(nargin < 2)
    m = 1;
end

srt = sqrt((1 + conj(x).*x./m));

f = sum(srt - 1);

% calculate gradient
g     = (x./srt)/m;

% hessian
h     = 1./(sqrt(m)*srt);

%
if isdistributed(f)
    f = gather(f);
end


