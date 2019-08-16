function A = R2mat(R,idx)
% convert  band storate format to sparse matrix, see also mat2R
%
% use:
%   A = R2mat(R,idx);
%   
% input:
%   R   - bands of the matrix
%   idx - indices such that R(i,j) = A(i,i+idx(j))
%
% output:
%   A - sparse matrix
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

m = size(R,1);
    
for k = 1:length(idx)
    i = idx(k);
    R(:,k) = [zeros(i,1);R(max(1-i,1):min(m-i,m),k);zeros(-i,1)];
end

A = spdiags(R,idx,m,m);