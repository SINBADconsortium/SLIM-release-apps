function [H] = HSS_marine(A, dim1, dim2, column, mode)

switch mode
   case 1  % FORWARD
      A = reshape(A, dim1, dim2);
  
      % creating the tensors needed to extract the area of interests.
      row = (opRestriction(dim1, column(1):column(2)));
      col = (opRestriction(dim2, column(3):column(4)));
      L = opKron(col,row);

      % now use the tensor on the data A to get the restricted domain that you want
      H = L*A(:);

      % reshape H
      length_row = length(column(1):column(2));  % length of the restricted data along the row
      length_col = length(column(3):column(4));  % length of the restricted data along the column
      H = reshape(H, length_row, length_col);
        
   otherwise
      A = reshape(A, length(column(1):column(2)), length(column(3):column(4)));
      row = (opRestriction(dim1, column(1):column(2)));
      col = (opRestriction(dim2, column(3):column(4)));
      L = opKron(col,row);
      H = L'*A(:);
      H = reshape(H, dim1, dim2);
end

end  % function

