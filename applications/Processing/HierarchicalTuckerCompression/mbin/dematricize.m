function X = dematricize(Xt,sz,row_dims,col_dims)
% Dematricize - restores the tensor shape to a matricized tensor
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   X = dematricize(Xt,sz,row_dims,{col_dims});
% 
% Input:
%   Xt       - input matrix of size sz(row_dims) x sz(col_dims)
%   sz       - size of the original tensor
%   row_dims - dimensions present along the rows of Xt
%   col_dims - dimensions present along the columns of Xt (default:
%              setdiff(1:length(sz),row_dims)
%   
% Output:
%   X        - dematricized tensor of size sz

if nargin < 3
    error('Requires at least three arguments');
end
if ~isfloat(Xt) || length(size(Xt)) ~=2
    error('Xt must be a real/complex matrix');
end
d = length(sz);
if exist('col_dims','var')==0
    col_dims = setdiff(1:length(sz),row_dims);    
else
    s = sort(union(row_dims,col_dims),'ascend');    
    assert(norm(s - (1:d)) == 0, 'row_dims and col_dims must be a partition of 1:d');
end

if prod(sz(row_dims)) ~= size(Xt,1) || prod(sz(col_dims)) ~= size(Xt,2)
    error('Incompatible row/column dimensions');
end

X = reshape(Xt,[sz(row_dims),sz(col_dims)]);

p([row_dims,col_dims]) = 1:d;

X = permute(X,p);

end