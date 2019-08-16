function y = partition(N,P,overlap)
% Partition a vector 1:N in to chunks of size P with a given overlap.
% Note, some combinations of N,P,overlap will not yield a partition that
% goes exactly up to N, so a row with a possibly different overlap could be added.
%
% Curt Da Silva, 2016
%
% Usage:
%   y = partition(N,P,overlap);
%
% Input:
%   N       - number of elements
%   P       - chunk size
%   overlap - number of elements to overlap
%  
% Output:
%   y       - number of partitions x P matrix of indices 
%             each row corresponds to a chunk of indices
%
    y = bsxfun(@plus,1:P,(0:(P-overlap):N-P)'); 
    if max(y(end,:))<N
        y = [y; (N-P+1):N];
    end
end