function [Q,R] = qr_std(A)
% QR_STD - Returns the (mathematically) standard QR components (R has
% strictly nonnegative diagonal) of the matrix A
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   [Q,R] = qr_std(A);
%
% Input:
%   A    - m x n matrix, m > n
%
% Output:
%   Q,R  - QR components of A
[Q,R] = qr(A,0);
D = diag(sign(diag(R)));
R = D * R; Q = Q * D;
end