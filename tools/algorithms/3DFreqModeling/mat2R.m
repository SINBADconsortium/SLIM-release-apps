function [R,idx] = mat2R(A)
% convert sparse matrix to band storate format
%
% use:
%   [R,idx] = mat2R(A);
%   
% input:
%   A - sparse matrix
%   
% output:
%   R   - bands of the matrix
%   idx - indices such that R(i,j) = A(i,i+idx(j))
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



%% get just the rows:
[R,idx] = spdiags(A);
for k = 1:length(idx)
    m = size(A,2);
    i = idx(k);
    R(:,k) = [zeros(-i,1);R(max(1+i,1):min(m+i,m),k);zeros(i,1)];
end