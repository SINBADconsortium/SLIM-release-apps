function B = padarray(A,pad)
%% PADARRAY - Pads a 2D or 3D array with zeros
% 
% Usage:
%   B = padarray(A,pad);   
%
% Input:
%   A   - 2d or 3d array to pad
%   pad - nonnegative integer vector of length == ndims(A), each
%         entry is the amount to pad in that dimension
% 
% Output:
%   B   - padded array
    
s = size(A);
if length(s) < length(pad)
   s = [s,ones(1,length(pad)-length(s))];
end
newdims = s + pad;

B = zeros(newdims);
switch length(newdims)
  case 2
    [I,J] = ndgrid(1:size(A,1),1:size(A,2));
    [I] = sub2ind(newdims,vec(I),vec(J));
  case 3
    [I,J,K] = ndgrid(1:size(A,1),1:size(A,2),1:size(A,3));
    I = sub2ind(newdims,vec(I),vec(J),vec(K));
end
B(I) = vec(A);

end
