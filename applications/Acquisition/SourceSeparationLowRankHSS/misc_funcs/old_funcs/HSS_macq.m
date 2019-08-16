function [H] = HSS_macq(A, dim1, dim2, column, mode)

%------------------------------------------------------------------
% HSS_macq performs HSS partitioning on the input matrix
%
% Use:
%   HSS_macq(A, dim1, dim2, column, mode)
%
% Input:   
%          A - input data matrix
%       dim1 - number of rows of input frequency slice
%       dim2 - number of columns of input frequency slice
%     column - from HSS_cell.m
%       mode - [1 : forward], [-1 : adjoint]

% Authors: Haneet Wason and Rajiv Kumar
%          Seismic Laboratory for Imaging and Modeling
%          Department of Earth, Ocean, and Atmospheric Sciences
%          The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------

switch mode
   case 1   % FORWARD
      A = reshape(A, dim1, dim2);
  
      % creating the tensors needed to extract the area of interests
      row = (opRestriction(dim1, column(1):column(2)));
      col = (opRestriction(dim2, column(3):column(4)));
      L = opKron(col,row);

      % use the tensor on the data A to get the restricted domain that you want
      H = L*A(:);

      % reshape H
      length_row = length(column(1):column(2));  % length of the restricted data along the row
      length_col = length(column(3):column(4));  % length of the restricted data along the column
      H = reshape(H, length_row, length_col);
        
   otherwise   % ADJOINT
      A = reshape(A, length(column(1):column(2)), length(column(3):column(4)));
      row = (opRestriction(dim1, column(1):column(2)));
      col = (opRestriction(dim2, column(3):column(4)));
      L = opKron(col,row);
      H = L'*A(:);
      H = reshape(H, dim1, dim2);
end

end  % function

