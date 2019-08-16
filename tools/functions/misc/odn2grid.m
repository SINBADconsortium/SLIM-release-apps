function varargout = odn2grid(o,d,n)
% produces grid from {o,d,n} info
% see also grid2odn, odnread, odnwrite
% 
% use:
%   [x,y,z] = odn2grid(o,d,n)
%
% input:
%   {o,d,n} - vectors describing grid in each dimension.
%
% output:
%   {x,y,...} - grid in each dimension: x = o(1) + [0:n(1)-1]*d(1), etc.
%

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


for k = 1:length(o)
    varargout{k} = o(k) + [0:n(k)-1]*d(k);
end
