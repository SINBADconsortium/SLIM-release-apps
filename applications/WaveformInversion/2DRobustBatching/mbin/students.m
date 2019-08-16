function [f,g,h] = students(x, df)
% students T misfit: .5*(1+df)*log(1+x.^2/df)
%
% use:
%   [f,g,h] = students(x, {df})
%
% input:
%   x   - vector
%   df  - degrees of freedom, default = 2
%
% output:
%   f   - student's T misfit
%   g   - gradient as vector with same size as x
%   h   - positive approximation of Hessian as vector with same size as x

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
 df = 2;
end

% calculate student's T gradient and positive definite Hessian
% approximation
f = 0.5 * sum((df + 1).*log(1+(conj(x(:)).*x(:))./df));
g = (x(:).*(df + 1))./(df + conj(x(:)).*x(:));
h = (df + 1)./(df + conj(x(:)).*x(:));

%
if isdistributed(f)
    f = gather(f);
end
