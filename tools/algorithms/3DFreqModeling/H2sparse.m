% Converts matrix from band-storage format to standard Matlab's sparse
% format.
%
% USE:
%   A = H2sparse(H,idx)
%
% INPUT:
%   H       - Matrix generated with helmholtz_3d_operto routine
%   idx     - Collection of indexes generated with helmholtz_3d_operto routine
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: April, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function A = H2sparse(H,idx)

   N = size(H,2);
   A = sparse(N,N);

   for k = 1:length(idx)
      i = idx(k);
      A = A + spdiags([zeros(i,1);transpose(H(k,max(1-i,1):min(N-i,N)));zeros(-i,1)],i,N,N);
   end

end
