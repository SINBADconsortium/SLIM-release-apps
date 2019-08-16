function y = ndind(x,totalNDim,dim,i)
%NDIND  nth-dimension indexing
%
%   ndind(x,totalNDim,dim,i) will return the i-th dimth-dimensional-slice of x
%
%   totalNDim   total number of dimensions in x
%   dim         desired dimension to slice
%   i           desired slice index in x

indBefore(1:(dim-1)) = {':'};
indAfter(1:(totalNDim-dim)) = {':'};
y = x(indBefore{:},i,indAfter{:});