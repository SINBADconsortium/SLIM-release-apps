%-------------------------------------------------------------------------------
% Estimates the source amplitude.
% 
% Setting hmin to -1 (or any negative number) will return the identity.
% 
% Use: 
%   w = src_estimate(Wh,Dpred,Dobs)
% Input:
%  Wh        - Weighting/mask matrix
%  Dpred     - One single slice of the predicted data
%  Dobs      - One singlie slice of the observed data
%
% Output:
%   w        - The scalar containing the source estimate correction.
% 
% Author: Tristan van Leween
%         Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: September, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-------------------------------------------------------------------------------
function w = src_estimate(Wh,Dpred,Dobs)
   w = ((Wh.*Dpred(:))'*(Wh.*Dobs(:)))/norm(Wh.*Dpred(:)).^2; 
   w = full(w);
end