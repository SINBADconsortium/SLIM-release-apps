function [U,B] = HOSVD(X,ranks)
% HOSVD - Higher Order SVD - quasi-optimal truncation of a tensor to Tucker
% format.
%  
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%
% Usage:
%   [U,B] = HOSVD(X,ranks);
% 
% Input:
%   X      - d-dimensional tensor to truncate
%   ranks  - positive integer vector of length == d of truncation ranks,
%            ranks(i) corresponds to dimension i
% Output:
%   [U,B]  - Tucker parameters
U = cell(length(ranks),1);
B = X;
for i=1:length(ranks)
   [W,~,~] = svd(matricize(X,i),'econ');
   U{i} = W(:,1:ranks(i));
   B = ttm(B,U{i}',i);
end

end