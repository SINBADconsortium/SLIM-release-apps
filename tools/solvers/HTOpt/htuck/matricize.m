function Xt = matricize(X,row_dims,col_dims)
%MATRICIZE - Matricizes a tensor along specified dimensions
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   Xt = matricize(X,row_dims,{col_dims})
%
% Input:
%   X        - d-dimensional array
%   row_dims - dimensions to place along the rows of the matricization
%   col_dims - dimensions to place along the columns of the matricization 
%              (default: setdiff(1:d,row_dims) )


if nargin < 2
    error('Need at least two arguments');
end

if exist('col_dims','var')==0
    col_dims = setdiff(1:length(size(X)),row_dims);
end


szX = size(X);
while length(szX) < max(union(row_dims,col_dims))
    szX = [szX, 1];
end
row_sz = prod(szX(row_dims)); col_sz = prod(szX(col_dims));

Xt = permute(X,[row_dims,col_dims]);
Xt = reshape(Xt,row_sz,col_sz);


end