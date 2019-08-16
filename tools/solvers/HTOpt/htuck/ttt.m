function Z = ttt(X,Y,x_dims,y_dims)
%TTT - Tensor tensor contraction
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%  Usage:
%    Z = ttt(X,Y)
%
%    Z = ttt(X,Y,x_dims,y_dims)
%
%  Input:
%    X - d-dimensional tensor of size n = n1 x ... x nd
%    Y - e-dimensional tensor of size m = m1 x ... x me
%  x_dims - dimensions (subset of 1:d) over which to contract X (default: 1:d) 
%  y_dims - dimensions (subset of 1:e) over which to contract Y (default: 1:e)
%  
%  Require:
%   length(x_dims) == length(y_dims) 
%     AND
%   size(X,x_dims(i)) == size(Y,y_dims(i)) for each i = 1, ..., length(x_dims)
%
%  Output:
%   Z - (d-length(x_dims))+(e-length(y_dims)) dimensional tensor of size 
%     n(x_dims^c) x m(y_dims^c) where m,n are as above and x_dims^c
%     = setdiff(1:d,x_dims) and likewise for y_dims^c
%
%  Note: If X is complex, will contract over conj(X)   
if nargin < 2
   error('Need at least two arguments'); 
end

if ~isfloat(X) || ~isfloat(Y)
    error('X,Y must be multidimensional real arrays');
end

d = length(size(X));
e = length(size(Y));

if exist('x_dims','var') == 0
    x_dims = 1:d;
end

if exist('y_dims','var') == 0
    y_dims = 1:e;
end

if length(x_dims) ~= length(y_dims)
    error('Number of contraction indices must be the same');
end

for i=1:length(x_dims)
    if size(X,x_dims(i)) ~= size(Y,y_dims(i))
        error('Tensor dimensions must agree over contraction indices');
    end
end

x_keepdims = setdiff(1:d,x_dims);
szX = size(X);
y_keepdims = setdiff(1:e,y_dims);
szY = size(Y);

Xmat = matricize(X,x_keepdims,x_dims);
Ymat = matricize(Y,y_keepdims,y_dims);

if ~isreal(Xmat)
    Z = conj(Xmat) * (Ymat.');
else
    Z = Xmat * Ymat';
end
Z = dematricize(Z,[szX(x_keepdims),szY(y_keepdims)],1:length(x_keepdims));


end