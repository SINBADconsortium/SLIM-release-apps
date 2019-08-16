function [f,g,h] = twonorms(x,s)
% least-squares misfit: .5*(x/s).^2
%
% use:
%   [f,g,h] = twonorms(x,{s})
%
% input:
%   x   - vector
%   s   - variance, either scalar or vector of same size as x, default = 1.
%
% output:
%   f   - least-squares misfit
%   g   - gradient as vector with same size as x
%   h   - Hessian as vector with same size as x

% Author: Tristan van Leeuwen
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
    s = 1;
end

f = .5*norm(x./s).^2;
g = x./s.^2;
h = ones(size(x))./s.^2;