function Ifull = idx1d2nd(dims,I)
% Converts a vectorized 1d index in to a n-dimensional index 
%
% Usage:
%   Ifull = idx1d2nd(dims,I);
%
% Input:
%   dims  - size of each dimension
%   I     - vectorized 1d index corresponding to dims
%
% Output:
%   Ifull - d x length(I) matrix, Ifull(:,j) = ind2sub(dims,I)
    d = length(dims);
    Ifull = zeros(length(I),d);
    Ifull_str = '[Ifull(:,1)';
    for i=2:d
       Ifull_str = [Ifull_str,',Ifull(:,' num2str(i) ')']; 
    end
    Ifull_str = [Ifull_str,'] ='];
    eval([Ifull_str ' ind2sub(dims,I);']);
    Ifull = Ifull';
end