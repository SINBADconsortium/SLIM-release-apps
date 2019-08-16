function [R,idx] = Rtransp(R,idx)
% transpose band-storage matrix
%
% use:
%   [R,idx] = Rtransp(R,idx)
%
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

n = size(R,1);

for k = 1:length(idx)  
    R(:,k) = circshift(R(:,k),idx(k));  
    idx(k) = -idx(k);   
end

