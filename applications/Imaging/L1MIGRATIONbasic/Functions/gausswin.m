function [w] = gausswin(M,alfa)
% this is function which can generate an gauss window
% Iuput: 
%	 M: Number of the window
%	 alfa: Exponent
% Author: Xiang Li



n = -(M-1)/2 : (M-1)/2;
w = exp((-1/2) * (alfa * n/((M-1)/2)) .^ 2)';