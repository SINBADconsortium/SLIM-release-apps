function map = seiscol(s,p)

%------------------------------------------------------------------------------
% seismic colormap [blue - red] 
%
% Use:
%   map = seiscol({s},{p})
%
% Input: 
%     s - center of colormap between [-1 1], default = 0 (center).
%     p - power < 1 to emphasize large values, > 1 to emphasize small values. default is 1.
% 
% Output:
%      map - colormap

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: September, 2011

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------------------

% number of values = 2*n + 1
n = 64;

% input checking
if (nargin < 1) 
    s = 0;
    p = 1;
elseif (nargin < 2)
    p = 1;
end

% construct colormap
a = linspace(0,1-1/(n-1),n*(1+s));
b = linspace(1-1/(n-1),0,n*(1-s));

map = [a 1 1+0*b; a 1 b; 1+0*a 1 b]'.^p;

end % function

