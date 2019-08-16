function R = bitrevorder(X)
% (internal function for FNT files)
% Purpose: rearrange vector X (or columns of arrays) to reverse bit order, upto
% max 2^k size <= length(X).
% Example :
%  Indices of a 4-length vector X (with first index at 0), namely     
%   k=0 -> 00 -> flip -> 00
%   k=1 -> 01 -> flip -> 10
%   k=2 -> 10 -> flip -> 01
%   k=3 -> 11 -> flip -> 11
% are reordered as 
%   k=0 -> 00 -> flip -> 00
%   k=2 -> 10 -> flip -> 01
%   k=1 -> 01 -> flip -> 10
%   k=3 -> 11 -> flip -> 11
% so that the reversed bit code of indices is now ordered.  
    
if (isvector(X))
    R = bitrevorder_vector(X);
else
    R = zeros(size(X));
    for k = 1:size(X,2)
        R(:,k) = bitrevorder_vector(X(:,k));
    end
end

function R = bitrevorder_vector(X)
[f,e] = log2(length(X));
I = dec2bin(0:pow2(0.5,e)-1);
R = X(bin2dec(I(:,e-1:-1:1)) + 1);
