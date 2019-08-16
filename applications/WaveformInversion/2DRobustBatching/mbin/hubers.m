function [f g h] = hubers(x, m)
% Huber misfit: abs(x) for x<=m; x.^2 for x>m
%
% use:
%   [f,g,h] = hubers(x,{m})
%
% input:
%   x - vector  
%   m - scale: misfit is L2 for |x|<=m and L1 for |x|>m, default = 1
%
% output:
%   f - huber misfit
%   g - gradient as vector with same size as x
%   h - diagonal of Hessian as vector with same size as x

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


% calculate Huber norm
Ia = find(abs(x) <= m);
Ib = find(abs(x) >  m);

f = 0;
if ~isempty(Ia)
    f = .5*x(Ia)'*x(Ia)./m;
end
if ~isempty(Ib)
    f = f + sum(abs(x(Ib)) - .5*m);
end

% calculate gradient
g     = x./m;
g(Ib) = sign(x(Ib));

% hessian
h     = ones(size(x))./m;
h(Ib) = 0;

%
if isdistributed(f)
    f = gather(f);
end


