%-----------------------------------------------------------------------------
% Transpose band-storage matrix.
% ATTENTION: This only transposes; it does not conjugate!
%
% use:
%   [H,idx] = Htransp(H,idx)
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
%-----------------------------------------------------------------------------
function [H,idx] = Htransp(H,idx)

% If matrix is short-fat
if size(H,1)<size(H,2)
   for k = 1:length(idx)  
      H(k,:) = circshift(H(k,:),[0 idx(k)]);  
      idx(k) = -idx(k);   
   end
   
else
% If matrix is tall-skinny
   for k = 1:length(idx)  
      H(:,k) = circshift(H(:,k),idx(k));  
      idx(k) = -idx(k);   
   end
end

