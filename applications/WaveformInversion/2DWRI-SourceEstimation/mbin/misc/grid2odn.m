function [o,d,n] = grid2odn(varargin)
% produces grid info for writing .odn files.
% see also odn2grid, odnwrite, odnread
%
% use:
%   [o,d,n] = grid2odn(x,y,...)
%
% input:
%   {x,y} - vectors describing regular grid in each dimension
%
% output:
%   o - [x(1)        y(1)        ... ]
%   d - [x(2) - x(1) y(2) - y(1) ... ]
%   n - [length(x)   length(y)   ... ]

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

o = [];
d = [];
n = [];

for k=1:length(varargin)
    x = varargin{k};
    o = [o x(1)];
    n = [n length(x)];
    if n(end)>1
        d = [d x(2)-x(1)];
    else
        d = [d 1];
    end
end