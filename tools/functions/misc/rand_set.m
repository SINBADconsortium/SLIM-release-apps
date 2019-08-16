function y = rand_set(n,nmax,sort_incr)
% Rand_set - Generate random integers in bins.
%
% Curt Da Silva, 2016
%
% Usage:
%   y = rand_set(n,nmax,sort_incr);
%
% Input:
%   n    - numbins x 1 vector of # of samples in each bin
%   nmax - numbins x 1 vector of maximum size (integer) of each bin (maximum entry is sum(nmax))
%   sort_incr - if true, sort in increasing order (default: true)
% 
% Output:
%   y    - random binned integers
% 
    if exist('sort_incr','var')==0, sort_incr = true; end
    
    nbins = length(nmax);
    y = [];
    n_tot = [0,cumsum(nmax)];
    for i=1:nbins
        t = randperm(nmax(i),n(i));
        y = [y,n_tot(i)+t];
    end
    if sort_incr
        y = sort(y,'ascend');
    end
end