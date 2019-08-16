function [ y ] = kaiser_window( x, r, b )
%KAISER_WINDOW Kaiser finite impulse response (FIR) window 
%
% Curt Da Silva, 2015
%
% Usage:
%   y = kaiser_window(x, r, b);
%
% Input:
%   x - input grid
%   r - half-window length
%   b - kaiser window parameter
% 
% Output:
%   y - kaiser window evaluated on the grid x

y = zeros(size(x));
I = find(abs(x)<= r);
y(I) = besseli(0, b*((1-(x(I)/r).^2).^(1/2)))/besseli(0,b);


end

